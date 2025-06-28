// supplementary functions for all-vs-all overlap detection with interval-hash-augmented minimizers
//
// author: Hajime Suzuki
// date: 2024/2/18-2025/4/2
//
#include <math.h>
#include <stdlib.h>
#include "kvec.h"
#include "kalloc.h"
#include "mmpriv.h"
#include "khash.h"
#include "ksort.h"

int mm_ra_del_full_intl(const mm_idx_t *mi, int max_ovh, int min_intl, int qlen, int n_regs, mm_reg1_t *regs) {
	int i, j;
	for (i = j = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		int rlen = mi->seq[r->rid].len;

		// keep non-internal and long half-internal matches
		int s = r->qs < max_ovh || r->rs < max_ovh, e = r->qe+max_ovh > qlen || r->re+max_ovh > rlen;
		int is_full = s && e, is_half = s || e;
		if (is_full || (is_half && r->qe-r->qs >= min_intl && r->re-r->rs >= min_intl)) regs[j++] = regs[i];
	}
	return j;
}

static inline int32_t get_for_qpos(int32_t qlen, const mm128_t *a)
{
	int32_t x = (int32_t)a->y;
	int32_t q_span = a->y>>32 & 0xff;
	if (a->x>>63)
		x = qlen - 1 - (x + 1 - q_span); // revert the position to the forward strand of query
	return x;
}

static int get_mini_idx(int qlen, const mm128_t *a, int32_t n, const mm128_t *mv)
{
	int32_t x, L = 0, R = n - 1;
	x = get_for_qpos(qlen, a);
	while (L <= R) { // binary search
		int32_t m = ((uint64_t)L + R) >> 1;
		int32_t y = (int32_t)mv[m].y >> 1;
		if (y < x) L = m + 1;
		else if (y > x) R = m - 1;
		else return m;
	}
	return -1;
}

void mm_ra_est_err(const mm_idx_t *mi, int min_cnt, int qlen, int n_regs, mm_reg1_t *regs, const mm128_t *a, int n, mm128_t *mv, mm_ivh_idx_t *ivh_idx) {
	// calc dv:f. the algorithm is the same as mm_est_err; this uses mv instead of mini_pos. this accumulates hits to upper 32bits of mv.y for cal_purge_cov
	int i;
	uint64_t sum_k = 0;
	float avg_k;

	if (n == 0) return;
	for (i = 0; i < n; ++i) {
		sum_k += mv[i].x & 0xff;
		mv[i].y = (uint32_t)mv[i].y;
	}
	avg_k = (float)sum_k / n;

	for (i = 0; i < n_regs; ++i) {
		mm_reg1_t *r = &regs[i];
		int32_t st, en, n_match, n_flt, n_tot, l_ref; //, raw = 0;
		r->div = -1.0f;
		if (r->cnt == 0) continue;
		r->aux = st = en = get_mini_idx(qlen, r->rev? &a[r->as + r->cnt - 1] : &a[r->as], n, mv);
		if (st < 0) {
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] logic inconsistency in mm_est_err(). Please contact the developer.\n");
			continue;
		}
		l_ref = mi->seq[r->rid].len;
		mm_ivh_comp_hits_pileup(min_cnt, r->rev, qlen, r->cnt, &a[r->as], n, mv, ivh_idx, &n_flt, &n_tot, &n_match);
		r->frac_flt = (double)n_flt / (n_flt + n_tot);
		r->frac_hit = (double)n_match / (n_flt + n_tot);

		if (r->qs > avg_k && r->rs > avg_k) ++n_tot;
		if (qlen - r->qs > avg_k && l_ref - r->re > avg_k) ++n_tot;
		r->div = n_match >= n_tot? 0.0f : (float)(1.0 - pow((double)n_match / n_tot, 1.0 / avg_k));
	}
}

int mm_ra_select_sub_indv(void *km, const mm_idx_t *mi, int max_ovh, float pri_ratio, int best_n, int qlen, int n_regs, mm_reg1_t *regs) {
	int i, j, k, n = n_regs;
	mm128_t *a;
	a = (mm128_t*)kmalloc(km, n * sizeof(mm128_t));

	if (best_n > 0) --best_n;
	for (i = 0; i < n; ++i) {
		const mm_reg1_t *r = &regs[i];
		int rlen = mi->seq[r->rid].len, qs, qe, rs, re, l, score;
		qs = r->qs<max_ovh? 0 : r->qs-max_ovh, qe = r->qe+max_ovh>qlen? qlen : r->qe+max_ovh;
		rs = r->rs<max_ovh? 0 : r->rs-max_ovh, re = r->re+max_ovh>rlen? rlen : r->re+max_ovh;
		l = qe-qs<re-rs? qe-qs : re-rs;
		score = (int)((float)(r->p? r->score : r->score0) * 10000.0 / (float)l + 0.4999);
		a[i].x = ((uint64_t)r->rid<<32) | score;
		a[i].y = i;
	}
	radix_sort_128x(a, a+n);

	for (i = j = 0; i < n; ++i) {
		if (i == n-1 || (a[i].x>>32) != (a[i+1].x>>32)) {
			int s = j+best_n<i? i-best_n : j;
			for (k = i; k >= s; --k) if ((float)(int)a[k].x < pri_ratio*(float)(int)a[i].x) break;
			for (; k >= j; --k) regs[a[k].y].cnt = 0;
			j = i+1;
		}
	}
	for (i = j = 0; i < n; ++i) {
		if (regs[i].cnt != 0) regs[j++] = regs[i];
		else if (regs[i].p) free(regs[i].p);
	}
	kfree(km, a);
	return j;
}

// for debugging
void mm_ra_print_seeds(const mm_idx_t *mi, int n_segs, const int *qlens, int64_t n_a, const mm128_t *a, const char *qname) {
	(void)n_segs;
	int i, prev_dir_rid = -1;
	for (i = 0; i < n_a; ++i) {
		// sorted by rpos
		int dir_rid = a[i].x>>32;
		if (dir_rid != prev_dir_rid) {
			int rid = dir_rid<<1>>1;
			fprintf(stderr, "#ref\t%s\t%d\n", mi->seq[rid].name, mi->seq[rid].len);
			fprintf(stderr, "#query\t%s\t%d\n", qname, qlens[0]);
			prev_dir_rid = dir_rid;
		}
		fprintf(stderr, "%s\t%d\t%c\t%s\t%d\n", mi->seq[prev_dir_rid<<1>>1].name, (int)a[i].x, "+-"[(int)(a[i].x>>63)], qname, (int)a[i].y);
	}
}
