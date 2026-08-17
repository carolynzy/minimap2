#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>
<<<<<<< HEAD
#include <string.h>
#include <ctype.h>
=======
>>>>>>> 3c28777e7e2dcc90f825de1b9f17a89cca7d4452
#define __STDC_LIMIT_MACROS
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "ksw2.h"
#include "kvec.h"
#include "khash.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct {
	int32_t st, en, cnt;
	int32_t score:30, strand:2;
} mm_idx_intv1_t;

typedef struct mm_idx_intv_s {
	int32_t n, m;
	mm_idx_intv1_t *a;
} mm_idx_intv_t;

typedef struct mm_idx_jjump_s {
	int32_t n, m;
	mm_idx_jjump1_t *a;
} mm_idx_jjump_t;

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	uint32_t i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	if (mi->B) {
		for (i = 0; i < 1U<<mi->b; ++i) {
			free(mi->B[i].p);
			free(mi->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)mi->B[i].h);
		}
	}
	if (mi->spsc) free(mi->spsc);
	if (mi->I) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->I[i].a);
		free(mi->I);
	}
	if (mi->J) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->J[i].a);
		free(mi->J);
	}
<<<<<<< HEAD
	if (mi->is_ref) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->is_ref[i].a);
		free(mi->is_ref);
	}
=======
>>>>>>> 3c28777e7e2dcc90f825de1b9f17a89cca7d4452
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n)
{
	int mask = (1<<mi->b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>mi->b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int64_t n = 0, n1 = 0;
	uint32_t i;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1U<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	for (i = 0; i < 1U<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				sum += kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k)&1) ++n1;
			}
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %ld (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf; total length: %ld\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), (long)n, 100.0*n1/n, (double)sum / n, (double)len / sum, (long)len);
}

int mm_idx_index_name(mm_idx_t *mi)
{
	khash_t(str) *h;
	uint32_t i;
	int has_dup = 0, absent;
	if (mi->h) return 0;
	h = kh_init(str);
	for (i = 0; i < mi->n_seq; ++i) {
		khint_t k;
		k = kh_put(str, h, mi->seq[i].name, &absent);
		if (absent) kh_val(h, k) = i;
		else has_dup = 1;
	}
	mi->h = h;
	if (has_dup && mm_verbose >= 2)
		fprintf(stderr, "[WARNING] some database sequences have identical sequence names\n");
	return has_dup;
}

int mm_idx_name2id(const mm_idx_t *mi, const char *name)
{
	khash_t(str) *h = (khash_t(str)*)mi->h;
	khint_t k;
	if (h == 0) return -2;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

int mm_idx_getseq_rev(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	const mm_idx_seq_t *s;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	s = &mi->seq[rid];
	if (en > s->len) en = s->len;
	st1 = s->offset + (s->len - en);
	en1 = s->offset + (s->len - st);
	for (i = st1; i < en1; ++i) {
		uint8_t c = mm_seq4_get(mi->S, i);
		seq[en1 - i - 1] = c < 4? 3 - c : c;
	}
	return en - st;
}

int mm_idx_getseq2(const mm_idx_t *mi, int is_rev, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	if (is_rev) return mm_idx_getseq_rev(mi, rid, st, en, seq);
	else return mm_idx_getseq(mi, rid, st, en, seq);
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	if (n == 0) return INT32_MAX;
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
			assert(absent && j == start_a + n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				kh_val(h, itr) = p->y;
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				seq->is_alt = 0;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0)
				mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag&MM_I_HPC, &s->a);
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	mm_idx_post(pl.mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads) // a simpler interface; deprecated
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	fp = mm_bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, 14, flag, 1<<18, n_threads, UINT64_MAX);
	mm_bseq_close(fp);
	return mi;
}

mm_idx_t *mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n, const char **seq, const char **name)
{
	uint64_t sum_len = 0;
	mm128_v a = {0,0,0};
	mm_idx_t *mi;
	khash_t(str) *h;
	int i, flag = 0;

	if (n <= 0) return 0;
	for (i = 0; i < n; ++i) // get the total length
		sum_len += strlen(seq[i]);
	if (is_hpc) flag |= MM_I_HPC;
	if (name == 0) flag |= MM_I_NO_NAME;
	if (bucket_bits < 0) bucket_bits = 14;
	mi = mm_idx_init(w, k, bucket_bits, flag);
	mi->n_seq = n;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n, sizeof(mm_idx_seq_t)); // ->seq is allocated from km
	mi->S = (uint32_t*)calloc((sum_len + 7) / 8, 4);
	mi->h = h = kh_init(str);
	for (i = 0, sum_len = 0; i < n; ++i) {
		const char *s = seq[i];
		mm_idx_seq_t *p = &mi->seq[i];
		uint32_t j;
		if (name && name[i]) {
			int absent;
			p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
			strcpy(p->name, name[i]);
			kh_put(str, h, p->name, &absent);
			assert(absent);
		}
		p->offset = sum_len;
		p->len = strlen(s);
		p->is_alt = 0;
		for (j = 0; j < p->len; ++j) {
			int c = seq_nt4_table[(uint8_t)s[j]];
			uint64_t o = sum_len + j;
			mm_seq4_set(mi->S, o, c);
		}
		sum_len += p->len;
		if (p->len > 0) {
			a.n = 0;
			mm_sketch(0, s, p->len, w, k, i, is_hpc, &a);
			mm_idx_add(mi, a.n, a.a);
		}
	}
	free(a.a);
	mm_idx_post(mi, 1);
	return mi;
}

/*************
 * index I/O *
 *************/

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5], i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		if (mi->seq[i].name) {
			uint8_t l = strlen(mi->seq[i].name);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->seq[i].name, 1, l, fp);
		} else {
			uint8_t l = 0;
			fwrite(&l, 1, 1, fp);
		}
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[5], i;
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 5, fp) != 5) return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		fread(&l, 1, 1, fp);
		if (l) {
			s->name = (char*)kmalloc(mi->km, l + 1);
			fread(s->name, 1, l, fp);
			s->name[l] = 0;
		}
		fread(&s->len, 4, 1, fp);
		s->offset = sum_len;
		s->is_alt = 0;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
		fread(mi->S, 4, (sum_len + 7) / 8, fp);
	}
	return mi;
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= 4) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	mm_idx_t *mi;
	if (r->is_idx) {
		mi = mm_idx_load(r->fp.idx);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	if (mi) {
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
		mi->index = r->n_parts++;
	}
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}

#include <ctype.h>
#include <zlib.h>
#include "ksort.h"
#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

int mm_idx_alt_read(mm_idx_t *mi, const char *fn)
{
	int n_alt = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return -1;
	ks = ks_init(fp);
	if (mi->h == 0) mm_idx_index_name(mi);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *p;
		int id;
		for (p = str.s; *p && !isspace(*p); ++p) { }
		*p = 0;
		id = mm_idx_name2id(mi, str.s);
		if (id >= 0) mi->seq[id].is_alt = 1, ++n_alt;
	}
	mi->n_alt = n_alt;
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] found %d ALT contigs\n", __func__, n_alt);
	return n_alt;
}

/***************
 * BED reading *
 ***************/

#define sort_key_bed(a) ((a).st)
KRADIX_SORT_INIT(bed, mm_idx_intv1_t, sort_key_bed, 4)

#define sort_key_end(a) ((a).en)
KRADIX_SORT_INIT(end, mm_idx_intv1_t, sort_key_end, 4)

static mm_idx_intv_t *mm_idx_bed_read_core(const mm_idx_t *mi, const char *fn, int read_junc, int min_sc)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	mm_idx_intv_t *I;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	I = CALLOC(mm_idx_intv_t, mi->n_seq);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		mm_idx_intv_t *r;
		mm_idx_intv1_t t = {-1,-1,-1,-1,0};
		char *p, *q, *bl, *bs;
		int32_t i, id = -1, n_blk = 0;
		for (p = q = str.s, i = 0;; ++p) {
			if (*p == 0 || *p == '\t') {
				int32_t c = *p;
				*p = 0;
				if (i == 0) { // chr
					id = mm_idx_name2id(mi, q);
					if (id < 0) break; // unknown name; TODO: throw a warning
				} else if (i == 1) { // start
					t.st = atol(q); // TODO: watch out integer overflow!
					if (t.st < 0) break;
				} else if (i == 2) { // end
					t.en = atol(q);
					if (t.en < 0) break;
				} else if (i == 4) { // BED score
					t.score = *q >= '0' && *q <= '9'? atol(q) : -1;
				} else if (i == 5) { // strand
					t.strand = *q == '+'? 1 : *q == '-'? -1 : 0;
				} else if (i == 9) {
					if (!isdigit(*q)) break;
					n_blk = atol(q);
				} else if (i == 10) {
					bl = q;
				} else if (i == 11) {
					bs = q;
					break;
				}
				if (c == 0) break;
				++i, q = p + 1;
			}
		}
		if (id < 0 || t.st < 0 || t.st >= t.en) continue; // contig ID not found, or other problems
		if (min_sc > 0 && t.score < min_sc) continue;
		r = &I[id];
		if (i >= 11 && read_junc) { // BED12
			int32_t st, sz, en;
			st = strtol(bs, &bs, 10); ++bs;
			sz = strtol(bl, &bl, 10); ++bl;
			en = t.st + st + sz;
			for (i = 1; i < n_blk; ++i) {
				mm_idx_intv1_t s = t;
				if (r->n == r->m) {
					r->m = r->m? r->m + (r->m>>1) : 16;
					r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
				}
				st = strtol(bs, &bs, 10); ++bs;
				sz = strtol(bl, &bl, 10); ++bl;
				s.st = en, s.en = t.st + st;
				en = t.st + st + sz;
				if (s.en > s.st) r->a[r->n++] = s;
			}
		} else {
			if (r->n == r->m) {
				r->m = r->m? r->m + (r->m>>1) : 16;
				r->a = (mm_idx_intv1_t*)realloc(r->a, sizeof(*r->a) * r->m);
			}
			r->a[r->n++] = t;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return I;
}

static mm_idx_intv_t *mm_idx_bed_read_merge(const mm_idx_t *mi, const char *fn, int read_junc, int min_sc)
{
	long n = 0, n0 = 0;
	int32_t i;
	mm_idx_intv_t *I;
	I = mm_idx_bed_read_core(mi, fn, read_junc, min_sc);
	if (I == 0) return 0;
	for (i = 0; i < mi->n_seq; ++i) {
		int32_t j, j0, k;
		mm_idx_intv_t *intv = &I[i];
		n0 += intv->n;
		radix_sort_bed(intv->a, intv->a + intv->n); // sort by st
		for (j = 1, j0 = 0; j <= intv->n; ++j) { // sort by st and then by end
			if (j == intv->n || intv->a[j].st != intv->a[j0].st) {
				radix_sort_end(intv->a + j0, intv->a + j);
				j0 = j;
			}
		}
		for (j = 1, j0 = 0, k = 0; j <= intv->n; ++j) { // merge intervals with the same (st, en)
			if (j == intv->n || intv->a[j].st != intv->a[j0].st || intv->a[j].en != intv->a[j0].en) {
				intv->a[k] = intv->a[j0];
				intv->a[k++].cnt = j - j0;
				j0 = j;
			}
		}
		intv->a = REALLOC(mm_idx_intv1_t, intv->a, k);
		intv->n = intv->m = k;
		n += k;
	}
	if (mm_verbose >= 3)
		fprintf(stderr, "[%s] read %ld introns, %ld of which are non-redundant\n", __func__, n0, n);
	return I;
}

int mm_idx_bed_read(mm_idx_t *mi, const char *fn, int read_junc)
{
	if (mi->h == 0) mm_idx_index_name(mi);
	mi->I = mm_idx_bed_read_merge(mi, fn, read_junc, -1);
	return 0;
}

int mm_idx_bed_junc(const mm_idx_t *mi, int32_t ctg, int32_t st, int32_t en, uint8_t *s)
{
	int32_t i, left, right;
	mm_idx_intv_t *r;
	memset(s, 0, en - st);
	if (mi->I == 0 || ctg < 0 || ctg >= mi->n_seq) return -1;
	r = &mi->I[ctg];
	left = 0, right = r->n;
	while (right > left) {
		int32_t mid = left + ((right - left) >> 1);
		if (r->a[mid].st >= st) right = mid;
		else left = mid + 1;
	}
	for (i = left; i < r->n; ++i) {
		if (st <= r->a[i].st && en >= r->a[i].en && r->a[i].strand != 0) {
			if (r->a[i].strand > 0) {
				s[r->a[i].st - st] |= 1, s[r->a[i].en - 1 - st] |= 2;
			} else {
				s[r->a[i].st - st] |= 8, s[r->a[i].en - 1 - st] |= 4;
			}
		}
	}
	return left;
}

/*********************************
 * Reading junctions for jumping *
 *********************************/

#define sort_key_jj(a) ((a).off)
KRADIX_SORT_INIT(jj, mm_idx_jjump1_t, sort_key_jj, 4)

#define sort_key_jj2(a) ((a).off2)
KRADIX_SORT_INIT(jj2, mm_idx_jjump1_t, sort_key_jj2, 4)

static void sort_jjump(mm_idx_jjump_t *jj2)
{
	int32_t j0, j, k;
	if (jj2 == 0 || jj2->n == 0) return;
	radix_sort_jj(jj2->a, jj2->a + jj2->n);
	for (j0 = 0, j = 1; j <= jj2->n; ++j) {
		if (j == jj2->n || jj2->a[j0].off != jj2->a[j].off) {
			radix_sort_jj2(jj2->a + j0, jj2->a + j);
			j0 = j;
		}
	}
	// the actual merge
	for (j0 = 0, j = 1, k = 0; j <= jj2->n; ++j) {
		if (j == jj2->n || jj2->a[j0].off != jj2->a[j].off || jj2->a[j0].off2 != jj2->a[j].off2) {
			int32_t t, cnt = 0;
			uint16_t flag = 0;
			for (t = j0; t < j; ++t) cnt += jj2->a[t].cnt, flag |= jj2->a[t].flag;
			jj2->a[k] = jj2->a[j0];
			jj2->a[k].cnt = cnt;
			jj2->a[k++].flag = flag;
			j0 = j;
		}
	}
	jj2->n = k;
	jj2->a = REALLOC(mm_idx_jjump1_t, jj2->a, k);
}

static mm_idx_jjump_t *mm_idx_bed2jjump(const mm_idx_t *mi, const mm_idx_intv_t *I, uint16_t flag)
{
	int32_t i;
	mm_idx_jjump_t *J;
	J = CALLOC(mm_idx_jjump_t, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i) {
		int32_t j, k;
		const mm_idx_intv_t *intv = &I[i];
		mm_idx_jjump_t *jj = &J[i];
		jj->n = intv->n * 2;
		jj->a = CALLOC(mm_idx_jjump1_t, jj->n);
		for (j = k = 0; j < intv->n; ++j) {
			jj->a[k].off = intv->a[j].st, jj->a[k].off2 = intv->a[j].en, jj->a[k].cnt = intv->a[j].cnt, jj->a[k].strand = intv->a[j].strand, jj->a[k++].flag = flag;
			jj->a[k].off = intv->a[j].en, jj->a[k].off2 = intv->a[j].st, jj->a[k].cnt = intv->a[j].cnt, jj->a[k].strand = intv->a[j].strand, jj->a[k++].flag = flag;
		}
		sort_jjump(jj);
	}
	return J;
}

static mm_idx_jjump_t *mm_idx_jjump_merge(const mm_idx_t *mi, const mm_idx_jjump_t *J0, const mm_idx_jjump_t *J1)
{
	int32_t i;
	mm_idx_jjump_t *J2;
	J2 = CALLOC(mm_idx_jjump_t, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i) {
		int32_t j, k;
		const mm_idx_jjump_t *jj0 = &J0[i], *jj1 = &J1[i];
		mm_idx_jjump_t *jj2 = &J2[i];
		jj2->n = jj0->n + jj1->n;
		jj2->a = CALLOC(mm_idx_jjump1_t, jj2->n);
		for (j = k = 0; j < jj0->n; ++j) jj2->a[k++] = jj0->a[j];
		for (j = 0; j < jj1->n; ++j) jj2->a[k++] = jj1->a[j];
		sort_jjump(jj2);
	}
	return J2;
}

int mm_idx_jjump_read(mm_idx_t *mi, const char *fn, int flag, int min_sc)
{
	int32_t i, j, n_anno = 0, n_misc = 0;
	mm_idx_intv_t *I;
	mm_idx_jjump_t *J;
	if (mi->h == 0) mm_idx_index_name(mi);
	I = mm_idx_bed_read_merge(mi, fn, 1, min_sc);
	J = mm_idx_bed2jjump(mi, I, flag);
	for (i = 0; i < mi->n_seq; ++i) free(I[i].a);
	free(I);
	if (mi->J) {
		mm_idx_jjump_t *J2;
		J2 = mm_idx_jjump_merge(mi, mi->J, J);
		for (i = 0; i < mi->n_seq; ++i) {
			free(mi->J[i].a); free(J[i].a);
		}
		free(mi->J); free(J);
		mi->J = J2;
	} else mi->J = J;
	for (i = 0; i < mi->n_seq; ++i) {
		for (j = 0; j < mi->J[i].n; ++j)
			if (mi->J[i].a[j].flag & MM_JUNC_ANNO) ++n_anno;
			else ++n_misc;
	}
	if (mm_verbose >= 3)
		fprintf(stderr, "[%s] there are %d annotated and %d other splice positions in the index\n", __func__, n_anno, n_misc);
	return 0;
}

static int32_t mm_idx_jump_get_core(int32_t n, const mm_idx_jjump1_t *a, int32_t x) // similar to mm_idx_find_intv()
{
	int32_t s = 0, e = n;
	if (n == 0) return -1;
	if (x < a[0].off) return -1;
	while (s < e) {
		int32_t mid = s + (e - s) / 2;
		if (x >= a[mid].off && (mid + 1 >= n || x < a[mid+1].off)) return mid;
		else if (x < a[mid].off) e = mid;
		else s = mid + 1;
	}
	assert(0);
}

const mm_idx_jjump1_t *mm_idx_jump_get(const mm_idx_t *db, int32_t cid, int32_t st, int32_t en, int32_t *n)
{
	mm_idx_jjump_t *s;
	int32_t l, r;
	*n = 0;
	if (cid >= db->n_seq || cid < 0 || db->J == 0) return 0;
	if (en < 0 || en > db->seq[cid].len) en = db->seq[cid].len;
	s = &db->J[cid];
	if (s->n == 0) return 0;
	l = mm_idx_jump_get_core(s->n, s->a, st);
	r = mm_idx_jump_get_core(s->n, s->a, en);
	*n = r - l;
	return &s->a[l + 1];
}

/****************
 * splice score *
 ****************/

typedef struct mm_idx_spsc_s {
	uint32_t n, m;
	uint64_t *a; // pos<<56 | score<<1 | acceptor
} mm_idx_spsc_t;

int32_t mm_idx_spsc_read2(mm_idx_t *idx, const char *fn, int32_t max_sc, float scale)
{
	gzFile fp;
	kstring_t str = {0,0,0};
	kstream_t *ks;
	int32_t dret, j;
	int64_t n_read = 0;

	fp = fn && strcmp(fn, "-") != 0? gzopen(fn, "rb") : gzdopen(0, "rb");
	if (fp == 0) return -1;
	if (idx->h == 0) mm_idx_index_name(idx);
	if (max_sc > 63) max_sc = 63;
	idx->spsc = Kcalloc(0, mm_idx_spsc_t, idx->n_seq * 2);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		mm_idx_spsc_t *s;
		char *p, *q, *name = 0;
		int32_t i, type = -1, strand = 0, cid = -1, score = -1;
		int64_t pos = -1;
		for (i = 0, p = q = str.s;; ++p) {
			if (*p == '\t' || *p == 0) {
				int c = *p;
				*p = 0;
				if (i == 0) {
					name = q;
				} else if (i == 1) {
					pos = atol(q);
				} else if (i == 2) {
					strand = *q == '+'? 1 : '-'? -1 : 0;
				} else if (i == 3) {
					type = *q == 'D'? 0 : *q == 'A'? 1 : -1;
				} else if (i == 4) {
					score = atoi(q);
					break;
				}
				if (c == 0) break;
				q = p + 1, ++i;
			}
		}
		if (i < 4) continue; // not enough fields
		if (scale > 0.0f && scale < 1.0f)
			score = score > 0.0f? (int)(score * scale + .499) : (int)(score * scale - .499);
		if (score > max_sc) score = max_sc;
		if (score < -max_sc) score = -max_sc;
		cid = mm_idx_name2id(idx, name);
		if (cid < 0 || type < 0 || strand == 0 || pos < 0) continue; // FIXME: give a warning!
		s = &idx->spsc[cid << 1 | (strand > 0? 0 : 1)];
		Kgrow(0, uint64_t, s->a, s->n, s->m);
		if (pos > 0 && pos < idx->seq[cid].len) { // ignore scores at the ends
			s->a[s->n++] = (uint64_t)pos << 8 | (score + KSW_SPSC_OFFSET) << 1 | type;
			++n_read;
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	for (j = 0; j < idx->n_seq * 2; ++j) {
		mm_idx_spsc_t *s = &idx->spsc[j];
		if (s->n > 0)
			radix_sort_64(s->a, s->a + s->n);
	}
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] read %ld splice scores\n", __func__, (long)n_read);
	return 0;
}

int32_t mm_idx_spsc_read(mm_idx_t *idx, const char *fn, int32_t max_sc)
{
	return mm_idx_spsc_read2(idx, fn, max_sc, 1.0f);
}

static int32_t mm_idx_find_intv(int32_t n, const uint64_t *a, int64_t x)
{
	int32_t s = 0, e = n;
	if (n == 0) return -1;
	if (x < a[0]>>8) return -1;
	while (s < e) {
		int32_t mid = s + (e - s) / 2;
		if (x >= a[mid]>>8 && (mid + 1 >= n || x < a[mid+1]>>8)) return mid;
		else if (x < a[mid]>>8) e = mid;
		else s = mid + 1;
	}
	assert(0);
}

int64_t mm_idx_spsc_get(const mm_idx_t *db, int32_t cid, int64_t st, int64_t en, int32_t rev, uint8_t *sc)
{
	const mm_idx_spsc_t *s;
	if (cid >= db->n_seq || cid < 0 || db->spsc == 0) return -1;
	if (en < 0 || en > db->seq[cid].len) en = db->seq[cid].len;
	memset(sc, 0xff, en - st);
	s = &db->spsc[cid << 1 | (!!rev)];
	if (s->n > 0) {
		int32_t j, l, r;
		l = mm_idx_find_intv(s->n, s->a, st);
		r = mm_idx_find_intv(s->n, s->a, en);
		for (j = l + 1; j <= r; ++j) {
			int64_t x = (s->a[j]>>8) - st;
			uint8_t score = s->a[j] & 0xff;
			assert(x <= en - st);
			if (x == en - st) continue;
			if (sc[x] == 0xff || sc[x] < score) sc[x] = score;
		}
	}
	return en - st;
}
<<<<<<< HEAD


/*************************************
 * IS-aware reference interval index *
 *************************************/

int mm_idx_is_ref_bed_read(mm_idx_t *mi, const char *fn)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t i;
	if (mi->h == 0) mm_idx_index_name(mi);
	fp = gzopen(fn, "r");
	if (fp == 0) return -1;
	mi->is_ref = CALLOC(mm_idx_intv_t, mi->n_seq);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *p = str.s, *tok[3];
		int n = 0, id;
		while (*p == ' ' || *p == '\t') ++p;
		if (*p == 0 || *p == '#') continue;
		for (tok[n++] = p; *p && n < 3; ++p) {
			if (*p == '\t' || *p == ' ') {
				*p++ = 0;
				while (*p == '\t' || *p == ' ') ++p;
				if (*p) tok[n++] = p;
			}
		}
		if (n >= 3 && (id = mm_idx_name2id(mi, tok[0])) >= 0) {
			mm_idx_intv_t *r = &mi->is_ref[id];
			mm_idx_intv1_t t;
			memset(&t, 0, sizeof(t));
			t.st = atol(tok[1]); t.en = atol(tok[2]); t.cnt = 1;
			t.score = -1; t.strand = 0;
			if (t.st < 0 || t.st >= t.en) continue;
			if (r->n == r->m) {
				r->m = r->m? r->m + (r->m>>1) : 16;
				r->a = REALLOC(mm_idx_intv1_t, r->a, r->m);
			}
			r->a[r->n++] = t;
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	/* merge intervals within 20 bp on the same contig */
	for (i = 0; i < (int32_t)mi->n_seq; ++i) {
		int32_t j, k;
		mm_idx_intv_t *intv = &mi->is_ref[i];
		if (intv->n < 2) continue;
		radix_sort_bed(intv->a, intv->a + intv->n);
		for (j = 1, k = 0; j < intv->n; ++j) {
			if (intv->a[j].st - intv->a[k].en <= 20) {
				if (intv->a[j].en > intv->a[k].en)
					intv->a[k].en = intv->a[j].en;
			} else {
				intv->a[++k] = intv->a[j];
			}
		}
		intv->n = intv->m = k + 1;
		intv->a = REALLOC(mm_idx_intv1_t, intv->a, intv->n);
	}
		if (mm_verbose >= 3) {
			int32_t total = 0;
			for (i = 0; i < (int32_t)mi->n_seq; ++i) total += mi->is_ref[i].n;
			fprintf(stderr, "[M::%s] loaded 0-based half-open reference IS intervals from %s\n", __func__, fn);
		fprintf(stderr, "[M::is_ref_summary] source=bed intervals=%d\n", total);
		for (i = 0; i < (int32_t)mi->n_seq; ++i) {
			int32_t j;
			for (j = 0; j < mi->is_ref[i].n; ++j)
				fprintf(stderr, "[M::is_ref_interval] source=bed ref=%s start=%d end=%d\n",
						mi->seq[i].name, mi->is_ref[i].a[j].st, mi->is_ref[i].a[j].en);
		}
	}
	return 0;
}

int mm_idx_is_ref_overlap(const mm_idx_t *mi, int32_t rid, int32_t st, int32_t en, int32_t window, int32_t *is_st, int32_t *is_en)
{
	const mm_idx_intv_t *r;
	int32_t left, right;
	if (mi->is_ref == 0 || rid < 0 || rid >= (int32_t)mi->n_seq || en <= st) return 0;
	r = &mi->is_ref[rid];
	left = 0; right = r->n;
	while (left < right) {
		int32_t mid = left + ((right - left) >> 1);
		if (r->a[mid].en + window <= st) left = mid + 1;
		else right = mid;
	}
	if (left < r->n && r->a[left].st - window < en) {
		if (is_st) *is_st = r->a[left].st;
		if (is_en) *is_en = r->a[left].en;
		return 1;
	}
	return 0;
}

int mm_idx_is_ref_count_overlap(const mm_idx_t *mi, int32_t rid, int32_t st, int32_t en, int32_t window)
{
	const mm_idx_intv_t *r;
	int32_t i, n = 0;
	if (mi->is_ref == 0 || rid < 0 || rid >= (int32_t)mi->n_seq || en <= st) return 0;
	r = &mi->is_ref[rid];
	for (i = 0; i < r->n; ++i) {
		int32_t nb_st = r->a[i].st - window, nb_en = r->a[i].en + window;
		if (nb_st < en && nb_en > st) {
			++n;
			if (mm_verbose >= 4)
				fprintf(stderr, "[M::is_ref_overlap] %s:%d-%d neighborhood=%d-%d query_ref_window=%d-%d\n",
						mi->seq[rid].name, r->a[i].st, r->a[i].en, nb_st, nb_en, st, en);
		}
	}
	return n;
}

int mm_idx_is_ref_count(const mm_idx_t *mi, int32_t rid)
{
	if (mi->is_ref == 0 || rid < 0 || rid >= (int32_t)mi->n_seq) return 0;
	return mi->is_ref[rid].n;
}

int mm_idx_is_ref_at(const mm_idx_t *mi, int32_t rid, int32_t i, int32_t *st, int32_t *en)
{
	if (mi->is_ref == 0 || rid < 0 || rid >= (int32_t)mi->n_seq) return 0;
	if (i < 0 || i >= mi->is_ref[rid].n) return 0;
	if (st) *st = mi->is_ref[rid].a[i].st;
	if (en) *en = mi->is_ref[rid].a[i].en;
	return 1;
}

int mm_idx_is_ref_find(mm_idx_t *mi, const mm_mapopt_t *opt)
{
	int32_t rid, s, strand, total = 0;
	if (opt->n_is_seq <= 0) return -1;
	if (mi->is_ref == 0) mi->is_ref = CALLOC(mm_idx_intv_t, mi->n_seq);
	for (rid = 0; rid < (int32_t)mi->n_seq; ++rid) {
		uint8_t *ref;
		int32_t rlen = mi->seq[rid].len, pos;
		ref = CALLOC(uint8_t, rlen);
		mm_idx_getseq(mi, rid, 0, rlen, ref);
		for (s = 0; s < opt->n_is_seq; ++s) {
			const mm_is_seq_t *is = &opt->is_seq[s];
			if (is->len < opt->is_min_ins_len || is->len > rlen) continue;
			for (strand = 0; strand < 2; ++strand) {
				int32_t seed_len = is->len < 15? is->len : 15;
				for (pos = 0; pos <= rlen - is->len; ++pos) {
					int32_t i, m = 0, n = 0;
					if (memcmp(ref + pos, is->seq[strand], seed_len) != 0) continue;
					for (i = 0; i < is->len; ++i) {
						uint8_t rc = ref[pos + i], ic = is->seq[strand][i];
						if (rc < 4 && ic < 4) ++n, m += rc == ic;
					}
					if (n > 0 && (float)m / n >= opt->is_min_id) {
						mm_idx_intv_t *r = &mi->is_ref[rid];
						mm_idx_intv1_t t;
						t.st = pos; t.en = pos + is->len; t.cnt = 1; t.score = m; t.strand = strand == 0? 1 : -1;
						if (r->n == r->m) {
							r->m = r->m? r->m + (r->m>>1) : 16;
							r->a = REALLOC(mm_idx_intv1_t, r->a, r->m);
						}
						r->a[r->n++] = t; ++total;
					}
				}
			}
		}
		free(ref);
	}
	for (rid = 0; rid < (int32_t)mi->n_seq; ++rid) {
		int32_t j, k;
		mm_idx_intv_t *r = &mi->is_ref[rid];
		if (r->n == 0) continue;
		radix_sort_bed(r->a, r->a + r->n);
		for (j = 1, k = 0; j < r->n; ++j) {
			if (r->a[j].st - r->a[k].en <= 20) {
				if (r->a[j].en > r->a[k].en) r->a[k].en = r->a[j].en;
			} else r->a[++k] = r->a[j];
		}
		r->n = r->m = k + 1;
		r->a = REALLOC(mm_idx_intv1_t, r->a, r->n);
	}
	if (mm_verbose >= 3) {
		fprintf(stderr, "[M::%s] detected %d reference IS interval(s) from IS sequence(s)\n", __func__, total);
		for (rid = 0; rid < (int32_t)mi->n_seq; ++rid) {
			int32_t j;
			for (j = 0; j < mi->is_ref[rid].n; ++j)
				fprintf(stderr, "[M::is_ref_interval] source=fasta ref=%s start=%d end=%d\n",
						mi->seq[rid].name, mi->is_ref[rid].a[j].st, mi->is_ref[rid].a[j].en);
		}
	}
	return total > 0? 0 : -1;
}

int mm_idx_is_ref_fill_seqs(const mm_idx_t *mi, mm_mapopt_t *opt)
{
	int32_t rid, j, old_n, add_n = 0, k;
	if (mi->is_ref == 0 || opt->n_is_seq > 0) return 0;
	for (rid = 0; rid < (int32_t)mi->n_seq; ++rid)
		add_n += mi->is_ref[rid].n;
	if (add_n == 0) return -1;
	old_n = opt->n_is_seq;
	opt->is_seq = (mm_is_seq_t*)realloc(opt->is_seq, (old_n + add_n) * sizeof(mm_is_seq_t));
	memset(&opt->is_seq[old_n], 0, add_n * sizeof(mm_is_seq_t));
	for (rid = 0, k = old_n; rid < (int32_t)mi->n_seq; ++rid) {
		const mm_idx_intv_t *r = &mi->is_ref[rid];
		for (j = 0; j < r->n; ++j, ++k) {
			int32_t len = r->a[j].en - r->a[j].st, x;
			mm_is_seq_t *s = &opt->is_seq[k];
			s->len = len;
			s->name = (char*)malloc(strlen(mi->seq[rid].name) + 64);
			sprintf(s->name, "%s:%d-%d", mi->seq[rid].name, r->a[j].st, r->a[j].en);
			s->seq[0] = (uint8_t*)malloc(len);
			s->seq[1] = (uint8_t*)malloc(len);
			mm_idx_getseq(mi, rid, r->a[j].st, r->a[j].en, s->seq[0]);
			for (x = 0; x < len; ++x) {
				int c = s->seq[0][len - 1 - x];
				s->seq[1][x] = c < 4? 3 - c : 4;
			}
		}
	}
	opt->n_is_seq = old_n + add_n;
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s] extracted %d IS sequence(s) from reference BED intervals\n", __func__, add_n);
	return 0;
}
=======
>>>>>>> 3c28777e7e2dcc90f825de1b9f17a89cca7d4452
