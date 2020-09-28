/*
 * Homogeneous field search.
 * Searches a CSV database of points for a collection of triangles
 * suitable for homogeneous Ingress fielding.
 *
 * HOMOGENOUS FIELDS
 *
 * An homogeneous field plan of degree n ("Hn") is a field plan
 * that an Ingress player can use to construct an n-layer deep field
 * through splitting and layering of nodes within a single
 * enclosing triangle.
 *
 * The densest construction possible in the game is an H6.
 *
 * An H1 field plan or H1 for short is a simple, triangular field.
 * An Ingress player stands at verticies of the field plan (each vertex is a
 * real-world location) and "throws" a link (edge) to another vertex.
 * When a new link forms a triangle, the game creates a "field" (or two)
 * and the player and the team receive points.
 *
 * An H2 is constructed from an an H1 as follows. Find a new vertex C
 * within the H1 field plan. C is the common vertex of three new H1
 * fields (triangles) sharing an edge with the original enclosing triangle.
 * The three new fields co-layer with the original field's single layer,
 * and the total construction is called an H2. Each point under an H2 is
 * "homogenously" covered by two fields.
 * An H2 plan always consists of 4 verticies and 3 interior H1 fields.
 *
 * For an H3 plan the process continues, finding new verticies
 * of each interior H1 such that what was previously an interior
 * H1 becomes an new interior H2 plan. Eventually an H3 is formed of
 * 7 vertices, 3 interior H2 fields and 9 interior H1 fields.
 *
 * In general, an Hn field contains 3**(n-1) H1 fields and
 * has (3**(n-1) - 1)/2 interior vertices (plus the 3 bounding
 * verticies).
 *
 *  Plan  #interior_verticies
 *  H2      1
 *  H3      4
 *  H4     13
 *  H5     40
 *  H6    121
 *
 * SEARCH ALGORITHM
 *
 * This program searchs a fixed collection of verticies
 * for a field plan of the specified depth.
 * It works by:
 *   - enumerating all possible triangles (counterclockwise)
 *   - ignoring triangles that don't have the minimum interior
 *     vertext requirement
 *   - recursively choosing spliiting verticies until an Hn
 *     field plan is complete; then the enclosing triangle is
 *     printed with its field plan, and its estimated area.
 *
 * A field is modeled as an array of triangles, ordered
 * outside-in. It turns out that the three immediately
 * interior sub-triangles of the i'th triangle will be at
 * array elements i*3+1, i*3+2, i*3+3.
 *
 * The search algorithm recurses on the end segment of the array,
 * incrementing i by one each descent, much the same was as
 * a permuting or combinatoric search might.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <err.h>
#include <limits.h>
#include <unistd.h>
#include <stdbool.h>
#include <assert.h>

#ifndef HAVE_STRLCAT
static size_t
strlcat(char *d, const char *s, size_t dsz)
{
	while (dsz && *d)
		dsz--, d++;
	while (dsz > 1 && *s)
		dsz--, *d++ = *s++;
	if (dsz)
		*d = '\0';
	return 0; /* Note: nonstandard because unused */
}
#endif

#ifndef NDEBUG
static int verbose_level = 0;
# define verbose(...)  do { if (verbose_level)   fprintf(stderr, __VA_ARGS__); } while (0)
# define verbose2(...) do { if (verbose_level>1) fprintf(stderr, __VA_ARGS__); } while (0)
#else
# define verbose_level 0
# define verbose(...)
# define verbose2(...)
#endif

/* Convert an angle in microdegrees into a strict string form */
static const char *
angle_str(long a /* µ° */) {
	static char buf[32];
	unsigned int len = 0;

	if (a < 0) {
		buf[len++] = '-';
		a = -a;
	}
	len += snprintf(buf + len, sizeof buf - len, "%ld", a / 1000000);
	a %= 1000000;
	if (a) {
		buf[len++] = '.';
		for (long b = 1000000; a && (b /= 10); ) {
			buf[len++] = '0' + (char)(a / b);
			a = a % b;
		}
	}
	buf[len] = '\0';
	return buf;
}

/*------------------------------------------------------------
 * 2D Vector
 * - vector components are integers for speed
 * - portal lat/lng coordinates are trivially projected into
 *   the plane (cylindrical proj) with µ° precision.
 *   (So this will be inaccurate for physically large distances)
 *
 * Note: Vectors are small structure values, and do not need
 * dynamic memory management. They can be passed and returned
 * from functions like other small types.
 */
typedef struct vector {
	long x;	/* lat: µ°  -90000000 … +90000000 */
	long y;	/* lng: µ° -180000000 …+180000000 */
	const char *name;
} vector;

#define NULL_VECTOR (struct vector){ 0, 0 }

static vector
vector_sub(const vector a, const vector b) {		/* a-b */
	return (vector){ a.x - b.x, a.y - b.y };
}

static long
vector_wedge(const vector a, const vector b) {		/* a^b */
	return a.x * b.y - a.y * b.x;
}

#ifndef NDEBUG
static const char *		/* representative form for debug */
vector_str(const vector a) {
	static char buf[1024];

	snprintf(buf, sizeof buf, "<%s,%s,",
		a.name ? a.name : "", angle_str(a.x));
	strlcat(buf, angle_str(a.y), sizeof buf);
	strlcat(buf, ">", sizeof buf);
	buf[sizeof buf - 1] = '\0';
	return buf;
}
#endif

static const char *				/* for DrawTools */
vector_json(const vector a) {
	static char buf[1024];

	snprintf(buf, sizeof buf, "{\"lat\":%s,\"lng\":", angle_str(a.x));
	strlcat(buf, angle_str(a.y), sizeof buf);
	strlcat(buf, "}", sizeof buf);
	return buf;
}

static double
tri_area(const vector tri[static 3])
{
	/* If we imagine two sides of the triangle
	 * are co-planar 3D vectors, then their cross product
	 * will be the area of the parallelogram they
	 * span. Halve that to get the triangle area.
	 * For counterclockwise triangles, the magnitude
	 * will be negative, so we negate the product. */

	vector ab = vector_sub(tri[1], tri[0]);
	vector bc = vector_sub(tri[2], tri[1]);
	return -0.5e-6 * vector_wedge(ab, bc);
}

/*------------------------------------------------------------
 * Vector set
 * A vector set is a growable sequence of vectors.
 */
typedef struct vset {
	unsigned int len;	/* length: nr of readable elements */
	unsigned int cap;	/* capacity: nr of writeable elements */
	vector el[];
} vset;

static void
vset_setcap(vset **v, unsigned int newcap)
{
	*v = realloc(*v, sizeof (vset) + newcap * sizeof (vector));
	if (!*v)
		err(1, "realloc");
	(*v)->cap = newcap;
}

static vset *
vset_new(void) {
	vset *ret = NULL;

	vset_setcap(&ret, 0); /* allocates! */
	ret->len = 0;
	return ret;
}

static void
vset_clear(vset **vp)
{
	(*vp)->len = 0;
}

static void
vset_free(vset *vs)
{
	free(vs);
}

#define VSET_ALLOCATE_INCREMENT 1

static _Bool
vset_append(vset **vp, const vector v) {
	if ((*vp)->len == (*vp)->cap)
		vset_setcap(vp, (*vp)->cap + VSET_ALLOCATE_INCREMENT);
	(*vp)->el[(*vp)->len++] = v;
	return true;
}

/*------------------------------------------------------------
 * CSV file reading
 * Comma-separated value files have lines of fields separated by
 * commas. Fields are optionally double quoted with quote-doubling.
 */

/* Returns a field, removing quotes */
static const char *
csv_string_field(char **sp) {
	char *s;
	static char buf[1024];
	unsigned int len = 0;
	_Bool inquote = false;

	for (s = *sp; *s && (inquote || *s != ','); s++) {
		if (*s != '"') {
			if (len < sizeof buf - 2)
				buf[len++] = *s;
		} else if (inquote && s[1] == '"') {
			if (len < sizeof buf - 2)
				buf[len++] = '"';
			s++;
		} else {
			inquote = !inquote;
		}
	}
	if (*s == ',') s++;
	*sp = s;
	buf[len] = '\0';
	return buf;
}

/* Reads a signed angle and any trailing comma. Returns true on success and updates *sp */
static _Bool
csv_angle_field(char **sp, long *angle_return /* µ° */) {
	char *s = *sp;
	long val = 0;
	unsigned int digits = 0;
	unsigned int intlen = -1;
	_Bool negative;

	if (*s == '-') {
		negative = true;
		++s;
	} else
		negative = false;

	for (;; s++) {
		if (*s == '.') {
			if (intlen != -1)
				return false;
			intlen = digits;
		} else if (*s < '0' || (*s - '0') > 9) {
			break;
		} else if (val > INT_MAX/10 ||
		           (val == INT_MAX/10 &&
			    (*s - '0') > INT_MAX%10)) {
			return false; /* overflow */
		} else {
			val = val * 10 + (*s - '0');
			digits++;
		}
	}
	if (!digits)
		return false;	/* no digits, e.g. "." or "-." */
	if (*s == ',')
		s++;		/* consume following comma */
	else if (*s)
		return false;	/* bad number (exponent notation?) */

	if (intlen == -1)
		intlen = digits;
	while (digits < intlen + 6) {
		val *= 10;
		digits++;
	}
	while (digits > intlen + 6) {
		val = val / 10;	/* round towards 0 */
		digits--;
	}

	*angle_return = negative ? -val : val;
	*sp = s;
	return true;
}

static vset *
read_csv(FILE *f, const char *filename)
{
	vset *result = vset_new();
	size_t linesz = 0;
	ssize_t n;
	char *line = NULL;
	unsigned int lineno = 0;

	while ((n = getline(&line, &linesz, f)) >= 0) {
		++lineno;
		/* skip heading line */
		if (lineno == 1 && strncmp(line, "Name,", 5) == 0)
			continue;

		/* Remove trailing newline */
		while (n && (line[n-1] == '\n' || line[n-1] == '\r'))
			line[--n] = '\0';

		/* "Tingalpa Creek Reserve",-27.506334,153.17907,"http://...",6 */
		/* Skip first field (name) */
		char *s = line;

		vector v;
		const char *name = csv_string_field(&s);

#if 0
		char buf[1024];
		snprintf(buf, sizeof buf, ":%u <%s>", lineno, name);
		v.name = strdup(buf);
#else
		v.name = strdup(name);
#endif

		if (!csv_angle_field(&s, &v.x) ||
		    !csv_angle_field(&s, &v.y))
			errx(1, "%s:%u: bad angle", filename, lineno);
		if (!vset_append(&result, v))
			err(1, "couldn't append vector");
	}
	return result;
}

/*------------------------------------------------------------
 * Vector functions
 */

/* Returns true iff ac is a left twist relative to ab,
 * i.e. counter-clockwise.
 * Note: returns false when ab=ac or ab=0 or ac=0 */
static _Bool
is_leftward(const vector ab, const vector ac)
{
	return vector_wedge(ab, ac) < 0;
}

/* Returns the elements p of v such that a->p is left of a->b.
 * Caller must free the returned vset. */
static vset *
vset_left_of(const vector a, const vector b, const vset *v)
{
	vset *result = vset_new();
	const vector ab = vector_sub(b, a);

	verbose2("vset_left_of(%.10s,%.10s,<%u>) = [", a.name, b.name, v->len);
	assert(ab.x || ab.y); /* same as assert(a != b) */

	for (unsigned i = 0; i < v->len; i++) {
		const vector ap = vector_sub(v->el[i], a);
		if (is_leftward(ab, ap)) {
			verbose2("%s%.10s", result->len?",":"", v->el[i].name);
			vset_append(&result, v->el[i]);
		}
	}
	verbose2("]\n");
	return result;
}

/*
 * Triangle generator.
 * State for generating all counter-clockwise triangles from
 * points in the given set,
 */
struct tri_generator {
	const vset *vs;
	vset *inside;
   /* private: */
	_Bool ji;
	unsigned int i, j;	/* 0 <= i < j < |vs| */
	unsigned int k;		/* 0 <= k < |ji?jileft:ijleft| */
	vset *ijleft;		/* points left of i->j */
	vset *jileft;		/* points left of j->i */
};

static void
tri_generator_init(struct tri_generator *gen, const vset *vs)
{
	gen->vs = vs;
	gen->i = 0;
	gen->j = 0; /* will increment early to j=1 */
	gen->k = 0;
	gen->ji = true;
	gen->ijleft = vset_new();
	gen->jileft = vset_new();
	gen->inside = NULL;
}

/* Generates the next unique triangle, and its
 * interior vertex set.
 *
 * Returns false when there are no more triangles.
 * The sets gen->ijleft and gen->inside will be NULL
 * (free) when this happens.
 *
 * Returns true for each new triangle. The triangle can
 * be read from the generator state structure:
 *     gen->inside           - interior points of ijk
 *
 * A caller may 'take' or free gen->inside at any time,
 * as long as they assign NULL to gen->inside to indicate
 * that this function no longer needs to manage it.
 */
static _Bool
tri_generator_next(struct tri_generator * restrict gen, vector ret[static 3]) {
	const vector *el = gen->vs->el;
	const unsigned int len = gen->vs->len;
	vset *kset = gen->ji ? gen->jileft : gen->ijleft;

	if (gen->i == len)
		goto done;

	while (gen->k == kset->len) {
		gen->k = 0;
		if (!gen->ji) {
			gen->ji = true;
			kset = gen->jileft;
			continue;
		}

		if (++gen->j == len - 1) {
			if (++gen->i == len - 2)
				goto done;
			gen->j = gen->i + 1;
		}

		/* partition [j+1..len] into left and right of ij */
		vset_clear(&gen->ijleft);
		vset_clear(&gen->jileft);
		vector ij = vector_sub(el[gen->j], el[gen->i]);
		for (unsigned k = gen->j + 1; k < len; k++) {
			vector ik = vector_sub(el[k], el[gen->i]);
			vset_append(is_leftward(ij, ik) ? &gen->ijleft : &gen->jileft, el[k]);
		}
		gen->ji = false;
		kset = gen->ijleft;
	}

	if (!gen->ji) {
		ret[0] = el[gen->i];
		ret[1] = el[gen->j];
	} else {
		ret[0] = el[gen->j];
		ret[1] = el[gen->i];
	}
	ret[2] = kset->el[gen->k];

	vset_free(gen->inside);
	vset *tmp = vset_left_of(ret[1], ret[2], kset);
	gen->inside = vset_left_of(ret[2], ret[0], tmp);
	vset_free(tmp);

	gen->k++;
	return true;

done:
	vset_free(gen->ijleft); gen->ijleft = NULL;
	vset_free(gen->jileft); gen->jileft = NULL;
	vset_free(gen->inside); gen->inside = NULL;
	gen->i = len;
	return false;
}

/*------------------------------------------------------------
 * Homogenous field searcing
 *
 * A potential H field is stored as an array of triangles ("tri"s).
 * The first tri is the bounding triangle of the Hm field.
 * The number of elements in the array is (3**m - 1)/2.
 * In the array, a triangle at position i bounds the three
 * sub-triangles at array positions 3i+1, 3i+2, 3i+3.
 */

struct hfield {
	unsigned int ntri;
	struct tri {
		vector v[3];		/* corners of triangle */
		vset *inner;		/* inner points of triangle
					   (may be NULL for H1 tris */
		unsigned int i;		/* index of chosen centre in .inner
					 * (not valid for H1 tris) */
	} tri[];
};

#define MAXDEPTH 6
#define MAXTRI   (((3*3*3*3*3*3)-1)/2)	/* (3⁶-1)/2 */

static unsigned int maxtri(unsigned int depth)
{
	unsigned int ret = 1;

	while (depth--)
		ret *= 3;
	return (ret - 1)/2;
}

/* The minimum number of interior verticies required for a tri at index i.
 * It will be the centroid and the sum of the three subtri's minima.
 * Computed early */
static unsigned char min_interior[MAXTRI];

static unsigned char
init_min_interior(unsigned int depth, unsigned int i)
{
	return min_interior[i] =
	         depth == 1
	         ? 0
	         : 1
	           + init_min_interior(depth - 1, i * 3 + 1)
	           + init_min_interior(depth - 1, i * 3 + 2)
	           + init_min_interior(depth - 1, i * 3 + 3);
}

static struct hfield *
hfield_new(unsigned int depth,
	const vector a, const vector b, const vector c,
	vset *inner /*given*/)
{
	unsigned i;
	unsigned int ntri = maxtri(depth);
	struct hfield *h;

	h = malloc(sizeof *h + ntri * sizeof h->tri[0]);
	if (!h)
		err(1, "malloc");
	h->ntri = ntri;
	h->tri[0].v[0] = a;
	h->tri[0].v[1] = b;
	h->tri[0].v[2] = c;
	h->tri[0].inner = inner;
	for (i = 1; i < maxtri(depth - 1); i++)
		h->tri[i].inner = vset_new();
	for (; i < ntri; i++)
		h->tri[i].inner = NULL;
	return h;
}

static void
hfield_free(struct hfield *h)
{
	if (h) {
		for (unsigned i = 0; i < h->ntri; i++)
			vset_free(h->tri[i].inner);
		free(h);
	}
}

#ifndef NDEBUG
__attribute__((sentinel))
static const char *
polygon_json(const vector *v1, ...)
{
	va_list ap;
	static char buf[2048];

	snprintf(buf, sizeof buf, "{\"type\":\"polygon\",\"latLngs\":[%s", v1 ? vector_json(*v1) : "");
	if (v1) {
	    va_start(ap, v1);
	    const vector *v;
	    while ((v = va_arg(ap, const vector *))) {
		strlcat(buf, ",", sizeof buf);
		strlcat(buf, vector_json(*v), sizeof buf);
	    }
	}
	strlcat(buf, "]}", sizeof buf);
	return buf;
}
#endif

static void
print_solution(const struct hfield *h)
{
	printf("Found solution, area: %f\n", tri_area(h->tri[0].v));
	printf("[");
	for (unsigned i = 0; i < h->ntri; i++) {
		if (i) printf(",\n ");
		printf("{\"type\":\"polygon\",\"latLngs\":[");
		for (unsigned j = 0; j < 3; j++) {
			if (j) printf(",");
			printf("%s", vector_json(h->tri[i].v[j]));
		}
		printf("],\"color\":\"#a24ac3\"}");
	}
	printf("]\n");
	fflush(stdout);
}

/*
 * Recursively search all the subtriangles of tri[i].
 * Assumes:
 *   tri[i].v[] is populated
 *   tri[i].inner is populated
 *
 * The general process is to try each inner vertex as
 * a splitting vertex; the three sub-tris' inners are
 * computed and assigned.
 * Recurses on i+1 so that all split choices are tried.
 */
static _Bool
rsearch(unsigned int const i, struct hfield *h)
{
	struct tri * const tri = h->tri;
	const vector a = tri[i].v[0];
	const vector b = tri[i].v[1];
	const vector c = tri[i].v[2];
	const vset *inner = tri[i].inner;

	verbose2("%*srsearch(%u,) tri[%u]={%.10s,%.10s,%.10s",i,"",i,i,
		a.name, b.name, c.name);
	if (verbose_level && inner) {
		verbose2(",inner=[");
		for (unsigned l = 0; l < inner->len; l++)
			verbose2("%s%.10s", l ? "," : "", inner->el[l].name);
		verbose2("]");
	}
	verbose2("} min_interior[%u]=%u\n", i, min_interior[i]);

	if (!min_interior[i]) {
		/* We have reached the first H1 subtri. Because all previous i
		 * must have been H2 (with non-zero min_interior), they would
		 * have set all their sub H1s up correctly. So tri[] is complete,
		 * and we have a solution. */
		return true; /* terminate the chain of rsearch() now */
	}

	/* Allocate space for three sub triangles */
	struct tri *mab = &tri[3*i+1];
	struct tri *mbc = &tri[3*i+2];
	struct tri *mca = &tri[3*i+3];
	mab->v[1] = a; mab->v[2] = b;
	mbc->v[1] = b; mbc->v[2] = c;
	mca->v[1] = c; mca->v[2] = a;
	for (unsigned mi = 0; mi < inner->len; mi++) {
		/* Try each interior point m as the splitting
		 * point for the triangle ABC. That is we form
		 * sub-triangles MAB, MBC, MCA */
		vector m = inner->el[mi];

		/* For an ABC that is H3+, divide its inner vset
		 * into the inner vsets of H2+ subtris MAB,MBC,MCA.
		 * We'll distribute the results to the appropriate
		 * subtris' inners. */
		if (min_interior[i] > 1) {
			vector ma = vector_sub(a, m);
			vector mb = vector_sub(b, m);
			vector mc = vector_sub(c, m);
			vset_clear(&mab->inner);
			vset_clear(&mbc->inner);
			vset_clear(&mca->inner);
			/* For each point p inside ABC */
			for (unsigned pi = 0; pi < inner->len; pi++) {
				if (pi == mi)
					continue; /* that isn't p=m */
				const vector p = inner->el[pi];

				vector mp = vector_sub(p, m);
				_Bool aleft = is_leftward(ma, mp);
				_Bool bleft = is_leftward(mb, mp);
				_Bool cleft = is_leftward(mc, mp);
				if (aleft && !bleft)
					vset_append(&mab->inner, p);
				else if (bleft && !cleft)
					vset_append(&mbc->inner, p);
				else
					vset_append(&mca->inner, p);
			}
			/* Abort if any of the sub-tris have impossibly
			 * few members */
			if (mab->inner->len < min_interior[3*i+1] ||
			    mbc->inner->len < min_interior[3*i+2] ||
			    mca->inner->len < min_interior[3*i+3])
			{
				verbose2("%*s distributed {%u,%u,%u} but need {%u+,%u+,%u+}\n",
					i, "", mab->inner->len, mbc->inner->len, mca->inner->len,
					min_interior[3*i+1],min_interior[3*i+2],min_interior[3*i+3]);
				continue;
			}
		}

		/* Complete the three interior triangles MAB, MBC, MCA */
		mab->v[0] = mbc->v[0] = mca->v[0] = m;

		/* Now we can recurse! (Usually to a sibling sometimes to a descendent tri) */
		if (rsearch(i + 1, h))
			return true; /* stop when solution found */
	}
	return false; /* no solutions */
}

int
main(int argc, char *argv[])
{
	unsigned int depth = 5;
	const char *filename = "Portal_Export.csv";
	int ch;
	int error = 0;

	while ((ch = getopt(argc, argv, "d:v")) != -1)
		switch (ch) {
		case 'd':
			if (sscanf(optarg, "%u", &depth) != 1) {
				warnx("bad integer for -d, '%s'", optarg);
				error = 1;
			}
			break;
		case 'v':
#ifndef NDEBUG
			verbose_level++;
#else
			warnx("verbose not supported (NDEBUG)");
#endif
			break;
		default:
			error = 1;
		}
	if (optind < argc)
		filename = argv[optind++];

	/* Check usage */
	if (error || optind < argc) {
		fprintf(stderr, "usage: %s [-v[v]] [-d depth] [file.csv]\n",
		    argv[0]);
		exit(2);
	}

	if (depth > MAXDEPTH)
		errx(1, "maximum depth is %u", MAXDEPTH);
	init_min_interior(depth, 0);

	/* Load in the portals */
	vset *vs;
	if (strcmp(filename, "-") == 0) {
		filename = "<stdin>";
		vs = read_csv(stdin, filename);
	} else {
		FILE *f = fopen(filename, "r");
		if (!f)
			err(1, "%s", filename);
		vs = read_csv(f, filename);
		fclose(f);
	}

	verbose("%s: %u portals\n", filename, vs->len);
	if (vs->len < 3 + min_interior[0]) {
		warnx("%s: insufficient portals for an H%u", filename, depth);
		return 0;
	}

	for (unsigned int i = 0; i < vs->len; i++)
		verbose2("[%u] = %s\n", i, vector_str(vs->el[i]));

	/* Iterate over all triangles in the set just loaded */
	struct tri_generator gen;
	tri_generator_init(&gen, vs);

	_Bool more;
#pragma omp parallel private(more)
	do {
		struct hfield *h = NULL;

#pragma omp critical
		{
			vector bound[3];
			more = tri_generator_next(&gen, bound);
			if (more) {
				if (verbose_level) {
					verbose("%u/%u/%u:\n", gen.i, gen.j, gen.k);
					verbose("  <%.16s,%.16s,%.16s>\n", bound[0].name, bound[1].name, bound[2].name);
					verbose2("  %s\n", polygon_json(&bound[0], &bound[1], &bound[2], NULL));
					verbose2("  inside: %u [", gen.inside->len);
					for (unsigned l = 0; l < gen.inside->len; l++)
						verbose2("%s%.16s", l ? "," : "", gen.inside->el[l].name);
					verbose2("]\n");
				}
				if (gen.inside->len >= min_interior[0]) {
					h = hfield_new(depth, bound[0], bound[1], bound[2], gen.inside);
					gen.inside = NULL; /* because now h 'owns' it */
				}
			}
		}

		if (h) {
			if (rsearch(0, h)) {
#pragma				omp critical
				print_solution(h);
			}
			hfield_free(h);
		}
	} while (more);

}
