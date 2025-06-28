// interval-hash-augmented minimizers
//
// author: Hajime Suzuki
// date: 2024/2/18-2025/4/2
// see also: mm_idx_patch_ivh in index.c
//
#include <stdlib.h>
#include "kalloc.h"
#include "mmpriv.h"
#include "ksort.h"

#define sort_key_pos(a) ((a).y<<32)	  // ignore upper 32bits
KRADIX_SORT_INIT(pos, mm128_t, sort_key_pos, 8)

#define beg_x(_p)		(((int32_t)(_p)->y>>1)-((_p)->x&0xff)+1)
#define beg_s(_p)		(((_p)->y>>1)-((_p)->q_span)+1)
#define end_x(_p)		((int32_t)(_p)->y>>1)
#define end_y(_y)		(*(_y)>>1)
#define end_y2(_y)		(*(_y)<<1>>2)

#define MARGIN_N (3)
static uint64_t hash[] = {
	0x58ea1ee2,	// 1
	0x41fc3e80, // 1.25
	0x3462e86b, // 1.5
	0x4cbf6848, // 1.75
	0x7bf817f0, // 2
	0x19b6c2ea, // 2.25
	0x69d22ca3, // 2.5
	0x5c49da04, // 2.75
	0x0ad06df1, // 3
	0x2161a558, // 3.25
	0x297f67ac, // 3.5
	0x32c2ea11, // 3.75
	0x2d5b49ac, // 4
	0x155f803c, // 4.25
	0x1584e4b5, // 4.5
	0x7431ccd0, // 4.75
	0x23faf39d, // 5
	0x1f1f17ac, // 5.25
	0x57064bd2, // 5.5
	0x0f00cf1c, // 5.75
	0x43390b8b, // 6
	0x36cee8a8, // 6.25
	0x173a7857, // 6.5
	0x1862821e, // 6.75
	0x4c669812, // 7
	0x7643748c, // 7.25
	0x4d550e1c, // 7.5
	0x7a1d81ba, // 7.75
	0x675497e1, // 8
	0x16ede062, // 8.25
	0x1b6d09a3, // 8.5
	0x2fe1504d, // 8.75
	0x2fa2328b, // 9
};
static uint64_t bnd[3] = {
	0,          // singleton
	0xfd2adec3, // query boundary
	0xba102f14, // target boundary
};

static void cal_unit_intv(int n, mm_ivh_iv_t *v, uint32_t wing) {
	int i, j, k, w1 = 2*wing;
	uint32_t min;
	for (i = j = 0, min = 1<<31; i < n-1; ++i) {
		if (v[i].iv == 0) {
			v[i].min_iv = min;
			min = 1<<31;
			j = i+1;
		} else if (v[i].iv < min) {
			int skip_first;
			min = v[i].iv;
			if (j+w1+1 <= i) ++j;
			skip_first = j+w1 == i;
			for (k = j+skip_first; k < i+1; ++k) if (min < v[k].min_iv) v[k].min_iv = min;
		} else {
			if (j+w1+1 <= i) {
				if (v[j].iv == min) {
					min = 1<<31;
					for (k = j+1; k < i+1; ++k) if (v[k].iv < min) min = v[k].iv;
				}
				++j;
			}
			v[i].min_iv = min;
		}
	}
	v[n-1].min_iv = min;
}

static void compute_hash(int n, mm_ivh_iv_t *v, uint32_t wing) {
	int i, j, b, e, w2 = wing;
	for (i = b = e = 0; i < n; ++i) {
		int s = 0, is_rev = v[i].is_rev;
		uint64_t q[32], d = v[i].min_iv>>2, rc = d>>1, h = 0;

		while (e < i+w2 && e < n && v[e].iv != 0) ++e;
		if (b + w2 < i) ++b;
		if (!is_rev) s = b + w2 - i;
		else s = i + w2 - e;
		if ((v[i].min_iv>>31) == 0) {
			for (j = 0; j < e-b; ++j) q[j] = d>0? ((v[b+j].iv)+rc) / d - 4 : 0;
			for (j = 0; j < e-b; ++j, ++s) {
				uint64_t *p = &q[is_rev? e-b-j-1 : j];
				h ^= hash[*p>31? 31 : *p]<<s;
			}
		} else h = bnd[v[i].min_iv&3];
		v[i].min_iv = h&0xffffff;    // min_iv field reused for hash

		if (v[i].iv == 0) b = e = i+1;
	}
}

static void cal_intv_and_compute_hash(const mm_idx_t *mi, int n, const uint64_t *y, uint32_t wing, uint32_t max_ivh_span, uint32_t k, uint64_t *vv) {
	int i;
	mm_ivh_iv_t *v = (mm_ivh_iv_t*)vv, *p;

	for (i = 0; i < n-1; ++i) {
		int64_t b = end_y(&y[i]), e = end_y(&y[i+1]) - k, dist = e-b;
		dist = dist<0? 0 : dist, dist = dist>max_ivh_span? 0 : dist;
		p = &v[i];
		p->is_rev = y[i]&1, p->iv = dist, p->min_iv = 1<<31, p->unused = p->aux = 0;
	}
	p = &v[n-1];
	p->is_rev = y[i]&1, p->iv = 0, p->min_iv = 1<<31, p->unused = p->aux = 0;
	cal_unit_intv(n, v, wing);
	if (mi->skip_bnd) {
		// suppress hash value to zero if it's too close to boundaries. applied only to the query side. this prevents false hits around boundaries
		for (i = 0; i < n; ++i) {
			uint32_t e = (uint32_t)y[i]>>1, qlen = mi->seq[y[i]>>32].len, mlen = MARGIN_N * v[i].min_iv;
			if (mlen > max_ivh_span/2) mlen = max_ivh_span/2;
			if (e < mlen || e + mlen > qlen) v[i].min_iv = 1<<31 | 2;
		}
	}
	compute_hash(n, v, wing);
}

int mm_ivh_compute_hash(const mm_idx_t *mi, int n, const uint64_t *y, uint32_t wing, uint32_t max_ivh_span, uint32_t k, uint64_t *vv) {
	int i, n_hash;
	if (n < 2 || wing == 0) return 1;
	cal_intv_and_compute_hash(mi, n, y, wing, max_ivh_span, k, vv);    // saves hash in upper 32bit of vv

	for (i = 0; i < n; ++i) vv[i] = vv[i]>>32<<32 | i;
	radix_sort_64(vv, vv+n);
	for (i = 0, n_hash = 0; i < n; ++i) if (i == n-1 || vv[i]>>32 != vv[i+1]>>32) ++n_hash;
	return n_hash;
}

int mm_ivh_flt_rep(int n, uint64_t *y, uint32_t rep_flt_span, uint32_t max_rep) {
	// filter out locally too-frequent kmers
	uint64_t *p, *s, *e;
	if (n < 2 || rep_flt_span == 0) return n;
	for (p = s = e = y; p < &y[n]; ++p) {
		while (s < p && end_y2(s) + rep_flt_span/2 <= end_y2(p)) ++s;
		while (e < &y[n] && end_y2(p) + rep_flt_span/2 > end_y2(e)) ++e;
		if (e-s >= max_rep) *p |= 1ULL<<63;
	}
	for (p = s = y; p < &y[n]; ++p) if (!(*p>>63)) *s++ = *p;
	return s-y;
}

static void update_hash_and_idx(size_t base, size_t n, mm128_t *mv, mm_ivh_idx_t *v, uint32_t wing) {
	size_t i, s, e;
	for (i = s = e = 0; i < n; ++i) {
		uint64_t idx = mv[i].y>>32, is_brk = v[i].iv == 0;
		while (i+wing > e && e < n && v[e].iv != 0) ++e;
		if (s+wing < i) ++s;

		// update hash
		mv[i].x ^= (uint64_t)v[i].mini_idx<<EMB_SIG_SHIFT;  // insert lower 24bits of pattern signature to upper 24bits of the original hash
		mv[i].y = ((base+i)<<32) | (uint32_t)mv[i].y;       // overwrite upper 32bits; points to index

		// build index
		v[i].fc = i-s, v[i].rc = e-i; // v[i].iv preserved
		v[i].mini_idx = idx, v[i].is_first = i == 0;
		if (is_brk) s = e = i+1;
	}
}

static void blank_idx(size_t base, size_t n, mm128_t *mv, mm_ivh_idx_t *v) {
	size_t i;
	for (i = 0; i < n; ++i) {
		uint64_t idx = mv[i].y>>32;
		mv[i].y = ((base+i)<<32) | (uint32_t)mv[i].y;
		v[i].fc = v[i].rc = v[i].iv = 0;
		v[i].mini_idx = idx, v[i].is_first = i == 0;
	}
}

mm_ivh_idx_t *mm_ivh_patch_sketch(void *km, int n, mm128_t *mv, int qlen, uint32_t wing, uint32_t max_ivh_span, uint32_t rep_flt_span, uint32_t max_rep, int skip_bnd) {
	size_t i, last_i, j;
	mm_ivh_iv_t *v;

	v = (mm_ivh_iv_t*)kmalloc(km, n * sizeof(mm_ivh_iv_t));
	if (n < 2) {
		memset(v, 0, n * sizeof(mm_ivh_iv_t));
		return (mm_ivh_idx_t*)v;
	}

	// sort by minimizers for building index; lower 32bits take priority over upper 24bits (where interval hash is embedded)
	for (i = 0; i < n; ++i) {
		uint64_t x = mv[i].x, hl = (uint32_t)(x>>8), hu = x>>40, l = x&0xff;
		mv[i].x = hl<<32 | hu<<8 | l;   // swap lower/upper bits
		mv[i].y |= i<<32;               // insert query position for use in radix_sort_pos
	}
	radix_sort_128x(mv, mv + n);  // sort by minimizers
	for (i = 0; i < n; ++i) {
		uint64_t x = mv[i].x, hl = (uint32_t)(x>>32), hu = x>>8, l = x&0xff;
		mv[i].x = hu<<40 | hl<<8 | l;   // swap back lower/upper bits
	}

	// extract interval and compute hash values from them
	for (i = last_i = 0; i < n; ++i) {
		if (i == n-1 || mv[i].x>>8 != mv[i+1].x>>8) {
			mm128_t *mv2 = &mv[last_i];
			mm_ivh_iv_t *v2 = &v[last_i], *p;
			int n2 = i+1-last_i;
			if (n2 >= 2) {
				radix_sort_pos(mv2, mv2+n2);
				for (j = 0; j < n2-1; ++j) {
					int64_t b = end_x(&mv2[j]), e = beg_x(&mv2[j+1]), dist = e-b;
					dist = dist<0? 0 : dist, dist = dist>max_ivh_span? 0 : dist;
					p = &v2[j];
					p->is_rev = mv2[j].y&1, p->iv = dist, p->min_iv = 1<<31, p->unused = p->aux = 0;
				}
				p = &v2[n2-1];
				p->is_rev = mv2[j].y&1, p->iv = 0, p->min_iv = 1<<31, p->unused = p->aux = 0;
				cal_unit_intv(n2, v2, wing);
				if (skip_bnd) {
					// skip hash augmentation if it's too close to boundaries. applied only to the query side. this prevents false hits around boundaries
					for (j = 0; j < n2; ++j) {
						uint32_t b = beg_x(&mv2[j]), e = end_x(&mv2[j]), mlen = MARGIN_N * v2[j].min_iv;
						if (mlen > max_ivh_span/2) mlen = max_ivh_span/2;
						if (b < mlen || e + mlen > qlen) v2[j].min_iv = 1<<31 | 1;
					}
				}
				compute_hash(n2, v2, wing);
				update_hash_and_idx(last_i, n2, mv2, (mm_ivh_idx_t*)v2, wing);
			} else blank_idx(last_i, n2, mv2, (mm_ivh_idx_t*)v2);
			last_i = i+1;
		}
	}

	// filter out locally too-frequent minimizers
	radix_sort_128x(mv, mv+n);  // sort again by minimizers
	for (i = last_i = 0; i < n; ++i) {
		if (i == n-1 || mv[i].x>>8 != mv[i+1].x>>8) {
			mm128_t *mv2 = &mv[last_i], *p, *s, *e;
			int n2 = i+1-last_i;
			if (n2 >= 2) {
				radix_sort_pos(mv2, mv2+n2);
				for (p = s = e = mv2; p < &mv2[n2]; ++p) {
					while (s < p && beg_x(s) + rep_flt_span/2 <= beg_x(p)) ++s;
					while (e < &mv2[n2] && end_x(p) + rep_flt_span/2 > end_x(e)) ++e;
					if (e-s >= max_rep) p->y |= 1ULL<<62;
				}
			}
			last_i = i+1;
		}
	}
	radix_sort_pos(mv, mv + n);  // sort back by pos
	return (mm_ivh_idx_t*)v;
}

static uint64_t weight(const mm_seed_t *a) {
	return (((1ULL<<63) / a->n)>>34);
}

static void seed_select(int n, mm_seed_t *a, uint32_t tie_rescue_w, uint32_t mid_occ) {
	int i;
	if (n == 0) return;

	// filter out too-frequent hits as in the original minimap2
	for (i = 0; i < n; ++i) if (a[i].n > mid_occ) a[i].flt = 1;

	// scan the number of hits (coverages) in a sliding window and recover a hit with the lowest occurrences,
	// if it does not make a tie with many in the window.
	if (n < tie_rescue_w) return;
	for (i = 0; i < n-tie_rescue_w+1; ++i) {
		int j, wt, max_wt = 0, max_i = -1, tie = 0;
		for (j = i; j < i+tie_rescue_w; ++j) {
			wt = weight(&a[j]);
			if (wt > max_wt) max_wt = wt, max_i = j, tie = 0;
			else if (wt == max_wt) ++tie;
		}
		if (tie < tie_rescue_w/2) a[max_i].flt = 0;
	}
}

mm_seed_t *mm_ivh_collect_matches(void *km, int *n_m, int qlen, const char *qname, int max_occ, int tie_rescue_w, const mm_idx_t *mi, const mm128_v *mv, int64_t *n_a)
{
	size_t i, j, k;
	mm_seed_t *m;
	m = mm_seed_collect_all(km, mi, mv, n_m);
	seed_select(*n_m, m, tie_rescue_w, max_occ);

	// merge flt flags
	for (i = j = k = 0, *n_a = 0; i < mv->n; ++i) {
		if (mv->a[i].y>>62 || m[j].y != mv->a[i].y) continue;
		if (m[j].flt) mv->a[i].y |= 1ULL<<62;
		else *n_a += m[j].n, m[k++] = m[j];
		++j;
	}
	assert(j == *n_m);
	*n_m = k;

	for (i = 0; i < *n_m; ++i) m[i].y = (uint32_t)m[i].y;
	return m;
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

void mm_ivh_comp_hits_pileup(int min_cnt, int rev, int qlen, int cnt, const mm128_t *a, int n, mm128_t *mv, mm_ivh_idx_t *idx, int *n_flt, int *n_tot, int *n_match) {
	// treat minimizers hit if it's embedded in another minimizer as its wing pattern and it hits
	int i, j, st, en;
	st = get_mini_idx(qlen, rev? &a[cnt - 1] : a, n, mv);

	// mark hits
	for (i = st, j = 0; i < n && j < cnt; ++i) {
		int32_t q;
		q = get_for_qpos(qlen, rev? &a[cnt - 1 - j] : &a[j]);
		if (q == (int32_t)mv[i].y>>1) mv[i].y |= 1ULL<<63, ++j;
	}
	en = i;

	// copy hits to the index array; idx is sorted by the original minimizers
	for (i = 0; i < n; ++i) idx[i].aux = mv[idx[i].mini_idx].y>>63;
	// mark embedded minimizers hit
	for (i = 0; i < n; ++i) {
		if (idx[i].aux == 0) continue;
		for (j = i-idx[i].fc; j < i+idx[i].rc+1; ++j)
			idx[j].aux |= 2;
	}

	// copy back hits
	for (i = 0; i < n; ++i) mv[idx[i].mini_idx].y |= (uint64_t)idx[i].aux>>1<<63;
	// pileup
	for (i = st, *n_flt = *n_match = *n_tot = 0; i < en; ++i) {
		if ((mv[i].y>>62)&1) {  // flt flag at bit 62
			*n_flt += 1;
			continue;
		}
		*n_tot += 1;
		if (mv[i].y>>63 == 0) continue;
		*n_match += 1;
		if (cnt >= min_cnt && ((mv[i].y>>32)&0x3fffffff) < 0x3fffffff) mv[i].y += 1ULL<<32;
	}

	for (i = 0; i < n; ++i)
		mv[i].y = mv[i].y<<1>>1, idx[i].aux = 0;
}
