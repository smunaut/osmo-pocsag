/* Quick & Dirty POCSAG experiment */

/* (C) 2012 by Sylvain Munaut <tnt@246tNt.com>
 * All Rights Reserved
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static int
hamming_weight(uint32_t data)
{
	int i, c=0;
	for (i=0; i<32 && data; i++) {
		c += data & 1;
		data >>= 1;
	}
	return c;
}

/* ------------------------------------------------------------------------ */

/*
 * the code used by POCSAG is a (n=31,k=21) BCH Code with dmin=5,
 * thus it could correct two bit errors in a 31-Bit codeword.
 * It is a systematic code.
 * The generator polynomial is: 
 *   g(x) = x^10+x^9+x^8+x^6+x^5+x^3+1
 * The parity check polynomial is: 
 *   h(x) = x^21+x^20+x^18+x^16+x^14+x^13+x^12+x^11+x^8+x^5+x^3+1
 *   g(x) * h(x) = x^n+1
 */

#define BCH_POLY 	0x769
#define BCH_N    	31
#define BCH_K    	21

static inline uint8_t
even_parity(uint32_t x) 
{
	x ^= x >> 16;
	x ^= x >> 8;
	x ^= x >> 4;
	x &= 0xf;
	return (0x6996 >> x) & 1;
}

static uint32_t
bch_syndrome(uint32_t data) 
{
	uint32_t shreg = data >> 1; /* throw away parity bit */
	uint32_t mask = 1L << (BCH_N-1), coeff = BCH_POLY << (BCH_K-1);
	int n = BCH_K;

	for(; n > 0; mask >>= 1, coeff >>= 1, n--)
		if (shreg & mask)
			shreg ^= coeff;

	if (even_parity(data))
		shreg |= (1 << (BCH_N - BCH_K));

	return shreg;
}

static uint32_t
bch_fix(uint32_t data)
{
	uint32_t t;
	int i, j;

	for (i=0; i<32; i++) {
		t = data ^ (1<<i);
		if (!bch_syndrome(t))
			return t;
	}

	for (i=0; i<32; i++) {
		for (j=i+1; j<32; j++) {
			t = data ^ ((1<<i) | (1<<j));
			if (!bch_syndrome(t))
				return t;
		}
	}

	return data;
}

/* ------------------------------------------------------------------------ */

#define POCSAG_TXT_LEN	256

struct pocsag_txt
{
	uint32_t bits;
	int nb;
	int nc;
	char buf[POCSAG_TXT_LEN];
};

static void
pocsag_txt_init(struct pocsag_txt *pt)
{
	pt->bits = 0;
	pt->nb = 0;
	pt->nc = 0;
	memset(pt->buf, 0x00, POCSAG_TXT_LEN);
}

static void
pocsag_txt_feed(struct pocsag_txt *pt, uint32_t d)
{
	/* Add 20 bits */
	pt->bits = (pt->bits << 20) | d;
	pt->nb += 20;

	/* Read as much char as possible */
	while (pt->nb >= 7) {
		char c = (pt->bits >> (pt->nb - 7)) & 0x7f;
		pt->buf[pt->nc++] =
			((c & 0x01) << 6) |
			((c & 0x02) << 4) |
			((c & 0x04) << 2) |
			 (c & 0x08)       |
			((c & 0x10) >> 2) |
			((c & 0x20) >> 4) |
			((c & 0x40) >> 6);
		pt->nb -= 7;
	}
}

/* ------------------------------------------------------------------------ */

#define POCSAG_SYNC     0x7cd215d8
#define POCSAG_IDLE     0x7a89c197

#define SYNC_MAX	10
#define SYNC_DEC	1
#define SYNC_INC	2

struct pocsag
{
	int b;
	int sync;

	uint32_t cw;
	int cw_bits;
	int cw_num;

	int txt_active;
	struct pocsag_txt txt;
};

static void
pocsag_init(struct pocsag *ps)
{
	memset(ps, 0x00, sizeof(struct pocsag));
	pocsag_txt_init(&ps->txt);
}

static void
pocsag_rx_bit(struct pocsag *ps, uint8_t bit)
{
	uint32_t cw;
	uint32_t syndrome;

	/* One more bit */
	ps->cw = (ps->cw << 1) | bit;
	ps->b++;

	/* Search for sync ? */
	if (!ps->sync) {
		if (ps->cw != POCSAG_SYNC)
			return;

		printf("SYNC @%d\n", ps->b-1);

		ps->sync = SYNC_MAX;

		ps->cw = 0;
		ps->cw_bits = 0;
		ps->cw_num = 0;

		return;
	}

	/* CW boundary */
	if (++ps->cw_bits != 32)
		return;
	
	ps->cw_bits = 0;

	/* Get / Check final CW */
	cw = ps->cw;

	syndrome = bch_syndrome(cw);

	if (syndrome) {
		cw = bch_fix(ps->cw);
		syndrome = bch_syndrome(cw);
	}

	/* SYNC tracking */
	if (syndrome) {
		ps->sync -= SYNC_DEC;
		if (ps->sync <= 0) {
			ps->sync = 0;
			printf(" <lost sync>\n");
			return;
		}
	} else {
		ps->sync += SYNC_INC;
		if (ps->sync > SYNC_MAX)
			ps->sync = SYNC_MAX;
	}

	/* Debug output */
	printf("%2d %08x (%08x) %d %2s",
		ps->cw_num, cw, ps->cw, hamming_weight(cw ^ ps->cw),
		syndrome == 0 ? "OK" : " "
	);

	if (ps->cw_num == -1) {
		if (cw == POCSAG_SYNC)
			printf(" => SYNC\n");
		ps->cw_num++;
		return;
	}

	if (cw == POCSAG_IDLE) {
		printf(" => IDLE");
		if (!syndrome) {
			ps->txt_active = 0;
		}
	} else if (cw == POCSAG_SYNC)
		printf(" => SYNC");
	else if (cw & (1<<31))
	{
		static uint32_t ld;
		char *tbl = "0123456789*U -)(";
		uint32_t d = (cw >> 11) & 0xfffff;
		int i;
		printf(" __ ");

		for (i=0; i<5; i++) {
			printf("%c", tbl[(d >> ((4-i)*4)) & 0xf]);
		}

#if 0
		if (cw_num & 1) {
			uint64_t fd = (ld << 20) | d;

			printf("  ");

			for (i=0; i<5; i++) {
				printf("%c", tbl[(d >> (i*5)) & 0xf]);
			}

		} else {
			ld = d;
		}
#endif
		if (ps->txt_active) {
			pocsag_txt_feed(&ps->txt, d);
			printf(" - %d %d TXT: |%s|", ps->txt.nb, ps->txt.nc, ps->txt.buf);
		}
	}
	else
	{
		printf(" => Addr: %d", (cw >> 11) & 3);

		if (!syndrome) {
			pocsag_txt_init(&ps->txt);
			ps->txt_active = 1;
		} else if (ps->txt_active) {
			uint32_t d = (cw >> 11) & 0xfffff;
			pocsag_txt_feed(&ps->txt, d);
		}
	}

	printf("\n");

	ps->cw_num ++;
	if (ps->cw_num == 16)
		ps->cw_num = -1;
}

int main(int argc, char *argv[])
{
	struct pocsag _ps, *ps=&_ps;
	FILE *f;

	if (argc != 2) {
		fprintf(stderr, "Usage: %s file.bits\n", argv[0]);
		return -1;
	}

	f = fopen(argv[1], "r");
	if (!f) {
		fprintf(stderr, "[!] Unable to open input file\n");
		return -1;
	}

	pocsag_init(ps);

	while (!feof(f))
	{
		uint8_t c;
		if ( fread(&c, 1, 1, f) != 1 )
			break;
		pocsag_rx_bit(ps, !c);
	}

	fclose(f);

	return 0;
}
