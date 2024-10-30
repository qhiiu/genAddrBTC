#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <cstdint>
#include <inttypes.h>
#include <math.h>
#include <emmintrin.h>

using namespace std;

#define P2PKH  0
#define P2SH   1
#define BECH32 2

//============== RIPEMD160.H - start ==============================================================================
namespace _ripemd160
{
	void inline Initialize(uint32_t *s)	{
		s[0] = 0x67452301ul;
		s[1] = 0xEFCDAB89ul;
		s[2] = 0x98BADCFEul;
		s[3] = 0x10325476ul; 
		s[4] = 0xC3D2E1F0ul;
	}

	inline uint32_t _rotl(uint32_t x, uint8_t r){	asm("roll %1,%0" : "+r"(x) : "c"(r));	return x;	}

	#define ROL(x,n) _rotl(x,n)
	#define f1(x, y, z) (x ^ y ^ z)
	#define f2(x, y, z) ((x & y) | (~x & z))
	#define f3(x, y, z) ((x | ~y) ^ z)
	#define f4(x, y, z) ((x & z) | (~z & y))
	#define f5(x, y, z) (x ^ (y | ~z))

	// hiiu_Round = Round nhưng bị trùng nên đổi thành hiiu_Round
	#define hiiu_Round(a,b,c,d,e,f,x,k,r) \
	a = ROL(a + f + x + k, r) + e; \
	c = ROL(c, 10);

	#define R11(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f1(b, c, d), x, 0, r)
	#define R21(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f2(b, c, d), x, 0x5A827999ul, r)
	#define R31(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f3(b, c, d), x, 0x6ED9EBA1ul, r)
	#define R41(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f4(b, c, d), x, 0x8F1BBCDCul, r)
	#define R51(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f5(b, c, d), x, 0xA953FD4Eul, r)
	#define R12(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f5(b, c, d), x, 0x50A28BE6ul, r)
	#define R22(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f4(b, c, d), x, 0x5C4DD124ul, r)
	#define R32(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f3(b, c, d), x, 0x6D703EF3ul, r)
	#define R42(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f2(b, c, d), x, 0x7A6D76E9ul, r)
	#define R52(a,b,c,d,e,x,r) hiiu_Round(a, b, c, d, e, f1(b, c, d), x, 0, r)

	/** Perform a RIPEMD-160 transformation, processing a 64-byte chunk. */
	void Transform(uint32_t *s, const unsigned char *chunk)
	{
		uint32_t a1 = s[0], b1 = s[1], c1 = s[2], d1 = s[3], e1 = s[4];
		uint32_t a2 = a1, b2 = b1, c2 = c1, d2 = d1, e2 = e1;
		uint32_t w[16];
		memcpy(w, chunk, 16 * sizeof(uint32_t));

		R11(a1, b1, c1, d1, e1, w[0], 11);
		R12(a2, b2, c2, d2, e2, w[5], 8);
		R11(e1, a1, b1, c1, d1, w[1], 14);
		R12(e2, a2, b2, c2, d2, w[14], 9);
		R11(d1, e1, a1, b1, c1, w[2], 15);
		R12(d2, e2, a2, b2, c2, w[7], 9);
		R11(c1, d1, e1, a1, b1, w[3], 12);
		R12(c2, d2, e2, a2, b2, w[0], 11);
		R11(b1, c1, d1, e1, a1, w[4], 5);
		R12(b2, c2, d2, e2, a2, w[9], 13);
		R11(a1, b1, c1, d1, e1, w[5], 8);
		R12(a2, b2, c2, d2, e2, w[2], 15);
		R11(e1, a1, b1, c1, d1, w[6], 7);
		R12(e2, a2, b2, c2, d2, w[11], 15);
		R11(d1, e1, a1, b1, c1, w[7], 9);
		R12(d2, e2, a2, b2, c2, w[4], 5);
		R11(c1, d1, e1, a1, b1, w[8], 11);
		R12(c2, d2, e2, a2, b2, w[13], 7);
		R11(b1, c1, d1, e1, a1, w[9], 13);
		R12(b2, c2, d2, e2, a2, w[6], 7);
		R11(a1, b1, c1, d1, e1, w[10], 14);
		R12(a2, b2, c2, d2, e2, w[15], 8);
		R11(e1, a1, b1, c1, d1, w[11], 15);
		R12(e2, a2, b2, c2, d2, w[8], 11);
		R11(d1, e1, a1, b1, c1, w[12], 6);
		R12(d2, e2, a2, b2, c2, w[1], 14);
		R11(c1, d1, e1, a1, b1, w[13], 7);
		R12(c2, d2, e2, a2, b2, w[10], 14);
		R11(b1, c1, d1, e1, a1, w[14], 9);
		R12(b2, c2, d2, e2, a2, w[3], 12);
		R11(a1, b1, c1, d1, e1, w[15], 8);
		R12(a2, b2, c2, d2, e2, w[12], 6);

		R21(e1, a1, b1, c1, d1, w[7], 7);
		R22(e2, a2, b2, c2, d2, w[6], 9);
		R21(d1, e1, a1, b1, c1, w[4], 6);
		R22(d2, e2, a2, b2, c2, w[11], 13);
		R21(c1, d1, e1, a1, b1, w[13], 8);
		R22(c2, d2, e2, a2, b2, w[3], 15);
		R21(b1, c1, d1, e1, a1, w[1], 13);
		R22(b2, c2, d2, e2, a2, w[7], 7);
		R21(a1, b1, c1, d1, e1, w[10], 11);
		R22(a2, b2, c2, d2, e2, w[0], 12);
		R21(e1, a1, b1, c1, d1, w[6], 9);
		R22(e2, a2, b2, c2, d2, w[13], 8);
		R21(d1, e1, a1, b1, c1, w[15], 7);
		R22(d2, e2, a2, b2, c2, w[5], 9);
		R21(c1, d1, e1, a1, b1, w[3], 15);
		R22(c2, d2, e2, a2, b2, w[10], 11);
		R21(b1, c1, d1, e1, a1, w[12], 7);
		R22(b2, c2, d2, e2, a2, w[14], 7);
		R21(a1, b1, c1, d1, e1, w[0], 12);
		R22(a2, b2, c2, d2, e2, w[15], 7);
		R21(e1, a1, b1, c1, d1, w[9], 15);
		R22(e2, a2, b2, c2, d2, w[8], 12);
		R21(d1, e1, a1, b1, c1, w[5], 9);
		R22(d2, e2, a2, b2, c2, w[12], 7);
		R21(c1, d1, e1, a1, b1, w[2], 11);
		R22(c2, d2, e2, a2, b2, w[4], 6);
		R21(b1, c1, d1, e1, a1, w[14], 7);
		R22(b2, c2, d2, e2, a2, w[9], 15);
		R21(a1, b1, c1, d1, e1, w[11], 13);
		R22(a2, b2, c2, d2, e2, w[1], 13);
		R21(e1, a1, b1, c1, d1, w[8], 12);
		R22(e2, a2, b2, c2, d2, w[2], 11);

		R31(d1, e1, a1, b1, c1, w[3], 11);
		R32(d2, e2, a2, b2, c2, w[15], 9);
		R31(c1, d1, e1, a1, b1, w[10], 13);
		R32(c2, d2, e2, a2, b2, w[5], 7);
		R31(b1, c1, d1, e1, a1, w[14], 6);
		R32(b2, c2, d2, e2, a2, w[1], 15);
		R31(a1, b1, c1, d1, e1, w[4], 7);
		R32(a2, b2, c2, d2, e2, w[3], 11);
		R31(e1, a1, b1, c1, d1, w[9], 14);
		R32(e2, a2, b2, c2, d2, w[7], 8);
		R31(d1, e1, a1, b1, c1, w[15], 9);
		R32(d2, e2, a2, b2, c2, w[14], 6);
		R31(c1, d1, e1, a1, b1, w[8], 13);
		R32(c2, d2, e2, a2, b2, w[6], 6);
		R31(b1, c1, d1, e1, a1, w[1], 15);
		R32(b2, c2, d2, e2, a2, w[9], 14);
		R31(a1, b1, c1, d1, e1, w[2], 14);
		R32(a2, b2, c2, d2, e2, w[11], 12);
		R31(e1, a1, b1, c1, d1, w[7], 8);
		R32(e2, a2, b2, c2, d2, w[8], 13);
		R31(d1, e1, a1, b1, c1, w[0], 13);
		R32(d2, e2, a2, b2, c2, w[12], 5);
		R31(c1, d1, e1, a1, b1, w[6], 6);
		R32(c2, d2, e2, a2, b2, w[2], 14);
		R31(b1, c1, d1, e1, a1, w[13], 5);
		R32(b2, c2, d2, e2, a2, w[10], 13);
		R31(a1, b1, c1, d1, e1, w[11], 12);
		R32(a2, b2, c2, d2, e2, w[0], 13);
		R31(e1, a1, b1, c1, d1, w[5], 7);
		R32(e2, a2, b2, c2, d2, w[4], 7);
		R31(d1, e1, a1, b1, c1, w[12], 5);
		R32(d2, e2, a2, b2, c2, w[13], 5);

		R41(c1, d1, e1, a1, b1, w[1], 11);
		R42(c2, d2, e2, a2, b2, w[8], 15);
		R41(b1, c1, d1, e1, a1, w[9], 12);
		R42(b2, c2, d2, e2, a2, w[6], 5);
		R41(a1, b1, c1, d1, e1, w[11], 14);
		R42(a2, b2, c2, d2, e2, w[4], 8);
		R41(e1, a1, b1, c1, d1, w[10], 15);
		R42(e2, a2, b2, c2, d2, w[1], 11);
		R41(d1, e1, a1, b1, c1, w[0], 14);
		R42(d2, e2, a2, b2, c2, w[3], 14);
		R41(c1, d1, e1, a1, b1, w[8], 15);
		R42(c2, d2, e2, a2, b2, w[11], 14);
		R41(b1, c1, d1, e1, a1, w[12], 9);
		R42(b2, c2, d2, e2, a2, w[15], 6);
		R41(a1, b1, c1, d1, e1, w[4], 8);
		R42(a2, b2, c2, d2, e2, w[0], 14);
		R41(e1, a1, b1, c1, d1, w[13], 9);
		R42(e2, a2, b2, c2, d2, w[5], 6);
		R41(d1, e1, a1, b1, c1, w[3], 14);
		R42(d2, e2, a2, b2, c2, w[12], 9);
		R41(c1, d1, e1, a1, b1, w[7], 5);
		R42(c2, d2, e2, a2, b2, w[2], 12);
		R41(b1, c1, d1, e1, a1, w[15], 6);
		R42(b2, c2, d2, e2, a2, w[13], 9);
		R41(a1, b1, c1, d1, e1, w[14], 8);
		R42(a2, b2, c2, d2, e2, w[9], 12);
		R41(e1, a1, b1, c1, d1, w[5], 6);
		R42(e2, a2, b2, c2, d2, w[7], 5);
		R41(d1, e1, a1, b1, c1, w[6], 5);
		R42(d2, e2, a2, b2, c2, w[10], 15);
		R41(c1, d1, e1, a1, b1, w[2], 12);
		R42(c2, d2, e2, a2, b2, w[14], 8);

		R51(b1, c1, d1, e1, a1, w[4], 9);
		R52(b2, c2, d2, e2, a2, w[12], 8);
		R51(a1, b1, c1, d1, e1, w[0], 15);
		R52(a2, b2, c2, d2, e2, w[15], 5);
		R51(e1, a1, b1, c1, d1, w[5], 5);
		R52(e2, a2, b2, c2, d2, w[10], 12);
		R51(d1, e1, a1, b1, c1, w[9], 11);
		R52(d2, e2, a2, b2, c2, w[4], 9);
		R51(c1, d1, e1, a1, b1, w[7], 6);
		R52(c2, d2, e2, a2, b2, w[1], 12);
		R51(b1, c1, d1, e1, a1, w[12], 8);
		R52(b2, c2, d2, e2, a2, w[5], 5);
		R51(a1, b1, c1, d1, e1, w[2], 13);
		R52(a2, b2, c2, d2, e2, w[8], 14);
		R51(e1, a1, b1, c1, d1, w[10], 12);
		R52(e2, a2, b2, c2, d2, w[7], 6);
		R51(d1, e1, a1, b1, c1, w[14], 5);
		R52(d2, e2, a2, b2, c2, w[6], 8);
		R51(c1, d1, e1, a1, b1, w[1], 12);
		R52(c2, d2, e2, a2, b2, w[2], 13);
		R51(b1, c1, d1, e1, a1, w[3], 13);
		R52(b2, c2, d2, e2, a2, w[13], 6);
		R51(a1, b1, c1, d1, e1, w[8], 14);
		R52(a2, b2, c2, d2, e2, w[14], 5);
		R51(e1, a1, b1, c1, d1, w[11], 11);
		R52(e2, a2, b2, c2, d2, w[0], 15);
		R51(d1, e1, a1, b1, c1, w[6], 8);
		R52(d2, e2, a2, b2, c2, w[3], 13);
		R51(c1, d1, e1, a1, b1, w[15], 5);
		R52(c2, d2, e2, a2, b2, w[9], 11);
		R51(b1, c1, d1, e1, a1, w[13], 6);
		R52(b2, c2, d2, e2, a2, w[11], 11);

		uint32_t t = s[0];
		s[0] = s[1] + c1 + d2;
		s[1] = s[2] + d1 + e2;
		s[2] = s[3] + e1 + a2;
		s[3] = s[4] + a1 + b2;
		s[4] = t + b1 + c2;
	}

} // namespace ripemd160

static const uint64_t sizedesc_32 = 32 << 3;
static const unsigned char pad[64] = { 0x80 };

void ripemd160_32(unsigned char *input, unsigned char *digest)
{
    uint32_t *s = (uint32_t *)digest;
    _ripemd160::Initialize(s);
    memcpy(input + 32, pad, 24);
    memcpy(input + 56, &sizedesc_32, 8);
    _ripemd160::Transform(s, input);
}

//============== RIPEMD160.H - end ===========================================================================
//============== SHA256 - start ==============================================================================
namespace _sha256
{
	static const unsigned char pad[64] = { 0x80 };

	#define _byteswap_ulong __builtin_bswap32
	#define _byteswap_uint64 __builtin_bswap64
	inline uint32_t _rotr(uint32_t x, uint8_t r){	asm("rorl %1,%0" : "+r"(x) : "c"(r));	return x;	}

	#define ROR(x,n) _rotr(x, n)
	#define S0(x) (ROR(x,2) ^ ROR(x,13) ^ ROR(x,22))
	#define S1(x) (ROR(x,6) ^ ROR(x,11) ^ ROR(x,25))
	#define s0(x) (ROR(x,7) ^ ROR(x,18) ^ (x >> 3))
	#define s1(x) (ROR(x,17) ^ ROR(x,19) ^ (x >> 10))

	#define Maj(x,y,z) ((x&y)^(x&z)^(y&z))

	// The following functions are equivalent to the above
	#define Ch(x,y,z) (z ^ (x & (y ^ z)))

	// SHA-256 round
	#define Round(a, b, c, d, e, f, g, h, k, w) \
		t1 = h + S1(e) + Ch(e,f,g) + k + (w); \
		t2 = S0(a) + Maj(a,b,c); \
		d += t1; \
		h = t1 + t2;

	#define WRITEBE32(ptr,x) *((uint32_t *)(ptr)) = _byteswap_ulong(x)
	#define WRITEBE64(ptr,x) *((uint64_t *)(ptr)) = _byteswap_uint64(x)
	#define READBE32(ptr) (uint32_t)_byteswap_ulong(*(uint32_t *)(ptr))

	// Initialise state
	void Initialize(uint32_t *s){
		s[0] = 0x6a09e667ul;
		s[1] = 0xbb67ae85ul;
		s[2] = 0x3c6ef372ul;
		s[3] = 0xa54ff53aul;
		s[4] = 0x510e527ful;
		s[5] = 0x9b05688cul;
		s[6] = 0x1f83d9abul;
		s[7] = 0x5be0cd19ul;
	}

	// Perform SHA-256 transformations, process 64-byte chunks
	void Transform(uint32_t *s, const unsigned char *chunk)
	{
		uint32_t t1, t2;
		uint32_t a = s[0], b = s[1], c = s[2], d = s[3], e = s[4], f = s[5], g = s[6], h = s[7];
		uint32_t w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15;

		Round(a, b, c, d, e, f, g, h, 0x428a2f98, w0 = READBE32(chunk + 0));
		Round(h, a, b, c, d, e, f, g, 0x71374491, w1 = READBE32(chunk + 4));
		Round(g, h, a, b, c, d, e, f, 0xb5c0fbcf, w2 = READBE32(chunk + 8));
		Round(f, g, h, a, b, c, d, e, 0xe9b5dba5, w3 = READBE32(chunk + 12));
		Round(e, f, g, h, a, b, c, d, 0x3956c25b, w4 = READBE32(chunk + 16));
		Round(d, e, f, g, h, a, b, c, 0x59f111f1, w5 = READBE32(chunk + 20));
		Round(c, d, e, f, g, h, a, b, 0x923f82a4, w6 = READBE32(chunk + 24));
		Round(b, c, d, e, f, g, h, a, 0xab1c5ed5, w7 = READBE32(chunk + 28));
		Round(a, b, c, d, e, f, g, h, 0xd807aa98, w8 = READBE32(chunk + 32));
		Round(h, a, b, c, d, e, f, g, 0x12835b01, w9 = READBE32(chunk + 36));
		Round(g, h, a, b, c, d, e, f, 0x243185be, w10 = READBE32(chunk + 40));
		Round(f, g, h, a, b, c, d, e, 0x550c7dc3, w11 = READBE32(chunk + 44));
		Round(e, f, g, h, a, b, c, d, 0x72be5d74, w12 = READBE32(chunk + 48));
		Round(d, e, f, g, h, a, b, c, 0x80deb1fe, w13 = READBE32(chunk + 52));
		Round(c, d, e, f, g, h, a, b, 0x9bdc06a7, w14 = READBE32(chunk + 56));
		Round(b, c, d, e, f, g, h, a, 0xc19bf174, w15 = READBE32(chunk + 60));

		Round(a, b, c, d, e, f, g, h, 0xe49b69c1, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0xefbe4786, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x0fc19dc6, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x240ca1cc, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x2de92c6f, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x4a7484aa, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x5cb0a9dc, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x76f988da, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0x983e5152, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0xa831c66d, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0xb00327c8, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0xbf597fc7, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0xc6e00bf3, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xd5a79147, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0x06ca6351, w14 += s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0x14292967, w15 += s1(w13) + w8 + s0(w0));

		Round(a, b, c, d, e, f, g, h, 0x27b70a85, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0x2e1b2138, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x4d2c6dfc, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x53380d13, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x650a7354, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x766a0abb, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x81c2c92e, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x92722c85, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0xa2bfe8a1, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0xa81a664b, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0xc24b8b70, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0xc76c51a3, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0xd192e819, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xd6990624, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0xf40e3585, w14 += s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0x106aa070, w15 += s1(w13) + w8 + s0(w0));

		Round(a, b, c, d, e, f, g, h, 0x19a4c116, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0x1e376c08, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x2748774c, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x34b0bcb5, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x391c0cb3, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x4ed8aa4a, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x5b9cca4f, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x682e6ff3, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0x748f82ee, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0x78a5636f, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0x84c87814, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0x8cc70208, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0x90befffa, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xa4506ceb, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0xbef9a3f7, w14 + s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0xc67178f2, w15 + s1(w13) + w8 + s0(w0));

		s[0] += a;
		s[1] += b;
		s[2] += c;
		s[3] += d;
		s[4] += e;
		s[5] += f;
		s[6] += g;
		s[7] += h;
	}

	void Transform2(uint32_t *s, const unsigned char *chunk)  // Compute SHA256(SHA256(chunk))[0]
	{
		uint32_t t1, t2;
		uint32_t w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15;

		uint32_t a = 0x6a09e667ul;
		uint32_t b = 0xbb67ae85ul;
		uint32_t c = 0x3c6ef372ul;
		uint32_t d = 0xa54ff53aul;
		uint32_t e = 0x510e527ful;
		uint32_t f = 0x9b05688cul;
		uint32_t g = 0x1f83d9abul;
		uint32_t h = 0x5be0cd19ul;

		Round(a, b, c, d, e, f, g, h, 0x428a2f98, w0 = READBE32(chunk + 0));
		Round(h, a, b, c, d, e, f, g, 0x71374491, w1 = READBE32(chunk + 4));
		Round(g, h, a, b, c, d, e, f, 0xb5c0fbcf, w2 = READBE32(chunk + 8));
		Round(f, g, h, a, b, c, d, e, 0xe9b5dba5, w3 = READBE32(chunk + 12));
		Round(e, f, g, h, a, b, c, d, 0x3956c25b, w4 = READBE32(chunk + 16));
		Round(d, e, f, g, h, a, b, c, 0x59f111f1, w5 = READBE32(chunk + 20));
		Round(c, d, e, f, g, h, a, b, 0x923f82a4, w6 = READBE32(chunk + 24));
		Round(b, c, d, e, f, g, h, a, 0xab1c5ed5, w7 = READBE32(chunk + 28));
		Round(a, b, c, d, e, f, g, h, 0xd807aa98, w8 = READBE32(chunk + 32));
		Round(h, a, b, c, d, e, f, g, 0x12835b01, w9 = READBE32(chunk + 36));
		Round(g, h, a, b, c, d, e, f, 0x243185be, w10 = READBE32(chunk + 40));
		Round(f, g, h, a, b, c, d, e, 0x550c7dc3, w11 = READBE32(chunk + 44));
		Round(e, f, g, h, a, b, c, d, 0x72be5d74, w12 = READBE32(chunk + 48));
		Round(d, e, f, g, h, a, b, c, 0x80deb1fe, w13 = READBE32(chunk + 52));
		Round(c, d, e, f, g, h, a, b, 0x9bdc06a7, w14 = READBE32(chunk + 56));
		Round(b, c, d, e, f, g, h, a, 0xc19bf174, w15 = READBE32(chunk + 60));

		Round(a, b, c, d, e, f, g, h, 0xe49b69c1, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0xefbe4786, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x0fc19dc6, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x240ca1cc, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x2de92c6f, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x4a7484aa, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x5cb0a9dc, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x76f988da, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0x983e5152, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0xa831c66d, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0xb00327c8, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0xbf597fc7, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0xc6e00bf3, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xd5a79147, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0x06ca6351, w14 += s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0x14292967, w15 += s1(w13) + w8 + s0(w0));

		Round(a, b, c, d, e, f, g, h, 0x27b70a85, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0x2e1b2138, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x4d2c6dfc, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x53380d13, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x650a7354, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x766a0abb, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x81c2c92e, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x92722c85, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0xa2bfe8a1, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0xa81a664b, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0xc24b8b70, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0xc76c51a3, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0xd192e819, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xd6990624, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0xf40e3585, w14 += s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0x106aa070, w15 += s1(w13) + w8 + s0(w0));

		Round(a, b, c, d, e, f, g, h, 0x19a4c116, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0x1e376c08, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x2748774c, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x34b0bcb5, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x391c0cb3, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x4ed8aa4a, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x5b9cca4f, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x682e6ff3, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0x748f82ee, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0x78a5636f, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0x84c87814, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0x8cc70208, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0x90befffa, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xa4506ceb, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0xbef9a3f7, w14 + s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0xc67178f2, w15 + s1(w13) + w8 + s0(w0));

		w0 = 0x6a09e667ul + a;
		w1 = 0xbb67ae85ul + b;
		w2 = 0x3c6ef372ul + c;
		w3 = 0xa54ff53aul + d;
		w4 = 0x510e527ful + e;
		w5 = 0x9b05688cul + f;
		w6 = 0x1f83d9abul + g;
		w7 = 0x5be0cd19ul + h;
		w8 = 0x80000000;
		w9 = 0;
		w10 = 0;
		w11 = 0;
		w12 = 0;
		w13 = 0;
		w14 = 0;
		w15 = 0x100;

		a = 0x6a09e667ul;
		b = 0xbb67ae85ul;
		c = 0x3c6ef372ul;
		d = 0xa54ff53aul;
		e = 0x510e527ful;
		f = 0x9b05688cul;
		g = 0x1f83d9abul;
		h = 0x5be0cd19ul;

		Round(a, b, c, d, e, f, g, h, 0x428a2f98, w0);
		Round(h, a, b, c, d, e, f, g, 0x71374491, w1);
		Round(g, h, a, b, c, d, e, f, 0xb5c0fbcf, w2);
		Round(f, g, h, a, b, c, d, e, 0xe9b5dba5, w3);
		Round(e, f, g, h, a, b, c, d, 0x3956c25b, w4);
		Round(d, e, f, g, h, a, b, c, 0x59f111f1, w5);
		Round(c, d, e, f, g, h, a, b, 0x923f82a4, w6);
		Round(b, c, d, e, f, g, h, a, 0xab1c5ed5, w7);
		Round(a, b, c, d, e, f, g, h, 0xd807aa98, w8);
		Round(h, a, b, c, d, e, f, g, 0x12835b01, w9);
		Round(g, h, a, b, c, d, e, f, 0x243185be, w10);
		Round(f, g, h, a, b, c, d, e, 0x550c7dc3, w11);
		Round(e, f, g, h, a, b, c, d, 0x72be5d74, w12);
		Round(d, e, f, g, h, a, b, c, 0x80deb1fe, w13);
		Round(c, d, e, f, g, h, a, b, 0x9bdc06a7, w14);
		Round(b, c, d, e, f, g, h, a, 0xc19bf174, w15);

		Round(a, b, c, d, e, f, g, h, 0xe49b69c1, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0xefbe4786, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x0fc19dc6, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x240ca1cc, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x2de92c6f, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x4a7484aa, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x5cb0a9dc, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x76f988da, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0x983e5152, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0xa831c66d, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0xb00327c8, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0xbf597fc7, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0xc6e00bf3, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xd5a79147, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0x06ca6351, w14 += s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0x14292967, w15 += s1(w13) + w8 + s0(w0));

		Round(a, b, c, d, e, f, g, h, 0x27b70a85, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0x2e1b2138, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x4d2c6dfc, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x53380d13, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x650a7354, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x766a0abb, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x81c2c92e, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x92722c85, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0xa2bfe8a1, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0xa81a664b, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0xc24b8b70, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0xc76c51a3, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0xd192e819, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xd6990624, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0xf40e3585, w14 += s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0x106aa070, w15 += s1(w13) + w8 + s0(w0));

		Round(a, b, c, d, e, f, g, h, 0x19a4c116, w0 += s1(w14) + w9 + s0(w1));
		Round(h, a, b, c, d, e, f, g, 0x1e376c08, w1 += s1(w15) + w10 + s0(w2));
		Round(g, h, a, b, c, d, e, f, 0x2748774c, w2 += s1(w0) + w11 + s0(w3));
		Round(f, g, h, a, b, c, d, e, 0x34b0bcb5, w3 += s1(w1) + w12 + s0(w4));
		Round(e, f, g, h, a, b, c, d, 0x391c0cb3, w4 += s1(w2) + w13 + s0(w5));
		Round(d, e, f, g, h, a, b, c, 0x4ed8aa4a, w5 += s1(w3) + w14 + s0(w6));
		Round(c, d, e, f, g, h, a, b, 0x5b9cca4f, w6 += s1(w4) + w15 + s0(w7));
		Round(b, c, d, e, f, g, h, a, 0x682e6ff3, w7 += s1(w5) + w0 + s0(w8));
		Round(a, b, c, d, e, f, g, h, 0x748f82ee, w8 += s1(w6) + w1 + s0(w9));
		Round(h, a, b, c, d, e, f, g, 0x78a5636f, w9 += s1(w7) + w2 + s0(w10));
		Round(g, h, a, b, c, d, e, f, 0x84c87814, w10 += s1(w8) + w3 + s0(w11));
		Round(f, g, h, a, b, c, d, e, 0x8cc70208, w11 += s1(w9) + w4 + s0(w12));
		Round(e, f, g, h, a, b, c, d, 0x90befffa, w12 += s1(w10) + w5 + s0(w13));
		Round(d, e, f, g, h, a, b, c, 0xa4506ceb, w13 += s1(w11) + w6 + s0(w14));
		Round(c, d, e, f, g, h, a, b, 0xbef9a3f7, w14 + s1(w12) + w7 + s0(w15));
		Round(b, c, d, e, f, g, h, a, 0xc67178f2, w15 + s1(w13) + w8 + s0(w0));

		s[0] = 0x6a09e667ul + a;
	}

} // namespace sha256

const uint8_t sizedesc_33[8] = { 0, 0, 0, 0, 0, 0, 1, 8 };
const uint8_t sizedesc_65[8] = { 0, 0, 0, 0, 0, 0, 2, 8 };

void sha256_33(unsigned char *input, unsigned char *digest) {
    uint32_t s[8];

    _sha256::Initialize(s);
    memcpy(input + 33, _sha256::pad, 23);
    memcpy(input + 56, sizedesc_33, 8);
    _sha256::Transform(s, input);

    WRITEBE32(digest, s[0]);
    WRITEBE32(digest + 4, s[1]);
    WRITEBE32(digest + 8, s[2]);
    WRITEBE32(digest + 12, s[3]);
    WRITEBE32(digest + 16, s[4]);
    WRITEBE32(digest + 20, s[5]);
    WRITEBE32(digest + 24, s[6]);
    WRITEBE32(digest + 28, s[7]);
}

void sha256_65(unsigned char *input, unsigned char *digest){
    uint32_t s[8];

    memcpy(input + 65, _sha256::pad, 55);
    memcpy(input + 120, sizedesc_65, 8);

    _sha256::Initialize(s);
    _sha256::Transform(s, input);
    _sha256::Transform(s, input + 64);

    WRITEBE32(digest, s[0]);
    WRITEBE32(digest + 4, s[1]);
    WRITEBE32(digest + 8, s[2]);
    WRITEBE32(digest + 12, s[3]);
    WRITEBE32(digest + 16, s[4]);
    WRITEBE32(digest + 20, s[5]);
    WRITEBE32(digest + 24, s[6]);
    WRITEBE32(digest + 28, s[7]);
}
 
void sha256_checksum(uint8_t *input, int length, uint8_t *checksum) {
    uint32_t s[8];
    uint8_t b[64];
    memcpy(b, input, length);
    memcpy(b + length, _sha256::pad, 56 - length);
    WRITEBE64(b + 56, length << 3);
    _sha256::Transform2(s, b);
    WRITEBE32(checksum, s[0]);

    // std::cout<<std::endl;  //solved problem 
}

//============== SHA256 - start ==============================================================================
//============== BASE58 - start ==============================================================================

/** All alphanumeric characters except for "0", "I", "O", and "l" */
static const char *pszBase58 = "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

std::string EncodeBase58(const unsigned char *pbegin, const unsigned char *pend)
{
    std::string ret;
    unsigned char digits[256];

    // Skip leading zeroes
    while (pbegin != pend && *pbegin == 0) {      ret.push_back('1');       pbegin++;     }
    int length = (int)(pend - pbegin);

    int digitslen = 1;
    digits[0] = 0;
    for (int i = 0; i < length; i++) {
        uint32_t carry = pbegin[i];
        for (int j = 0; j < digitslen; j++) {
            carry += (uint32_t)(digits[j]) << 8;
            digits[j] = (unsigned char)(carry % 58);
            carry /= 58;
        }
        while (carry > 0) {     digits[digitslen++] = (unsigned char)(carry % 58);     carry /= 58;     }
    }

    // reverse
    for (int i = 0; i < digitslen; i++){	ret.push_back(pszBase58[digits[digitslen - 1 - i]]);	}

    return ret;
}
//============== BASE58 - end ================================================================================
//============== INT.H - start ================================================================================
//---------- INT.H
#define NB64BLOCK 5
#define NB32BLOCK 10

class Int {
	public:
		Int();
		Int(int64_t i64);
		Int(uint64_t u64);
		Int(Int* a);

		void Add(uint64_t a);
		void Add(Int* a);
		void Add(Int* a, Int* b);
		void AddOne();
		void Sub(uint64_t a);
		void Sub(Int* a);
		void Sub(Int* a, Int* b);
		void Mult(Int* a);
		uint64_t Mult(uint64_t a);
		uint64_t Mult(Int* a, uint64_t b);
		uint64_t IMult(Int* a, int64_t b);
		void Mult(Int* a, Int* b);
		void Neg();

		// Comp 
		bool IsGreaterOrEqual(Int* a);
		bool IsEqual(Int* a);
		bool IsZero();
		bool IsOne();
		bool IsPositive();
		bool IsNegative();
		bool IsEven();

		// Setup field
		static void SetupField(Int* n, Int* R = NULL, Int* R2 = NULL, Int* R3 = NULL, Int* R4 = NULL);

		void ModInv();                             // this <- this^-1 (mod n)
		void MontgomeryMult(Int* a, Int* b);        // this <- a*b*R^-1 (mod n)
		void ModAdd(Int* a);                       // this <- this+a (mod n) [0<a<P]
		void ModAdd(Int* a, Int* b);                // this <- a+b (mod n) [0<a,b<P]
		void ModSub(Int* a);                       // this <- this-a (mod n) [0<a<P]
		void ModSub(Int* a, Int* b);               // this <- a-b (mod n) [0<a,b<P]
		void ModMul(Int* a, Int* b);                // this <- a*b (mod n) 
		void ModNeg();                             // this <- -this (mod n)

		// Specific SecpK1
		void ModMulK1(Int* a, Int* b);
		void ModMulK1(Int* a);
		void ModSquareK1(Int* a);

		// Size
		int GetSize();       // Number of significant 32bit limbs

		// Setter
		void SetInt32(uint32_t value);
		void Set(Int* a);
		void SetBase10(const char* value);
		void SetBase16(const char* value);
		void SetBaseN(int n, const char* charset, const char* value);

		// Getter
		unsigned char GetByte(int n);
		void Get32Bytes(unsigned char* buff);

		// To String
		std::string GetBase10();
		std::string GetBase16();
		std::string GetBaseN(int n, char* charset);

		union {
			uint32_t bits[NB32BLOCK];
			uint64_t bits64[NB64BLOCK];
		};

	private:
		void MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22, uint64_t* cu, uint64_t* cv);
		void MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22);
		uint64_t AddCh(Int* a, uint64_t ca, Int* b, uint64_t cb);
		uint64_t AddCh(Int* a, uint64_t ca);
		uint64_t AddC(Int* a);
		void AddAndShift(Int* a, Int* b, uint64_t cH);
		void CLEAR();
		void DivStep62(Int* u, Int* v, int64_t* eta, int* pos, int64_t* uu, int64_t* uv, int64_t* vu, int64_t* vv);
};

// Missing intrinsics
static uint64_t inline _umul128(uint64_t a, uint64_t b, uint64_t* h) {
	uint64_t rhi;
	uint64_t rlo;
	__asm__("mulq  %[b];" :"=d"(rhi), "=a"(rlo) : "1"(a), [b]"rm"(b));
	*h = rhi;
	return rlo;
}

static int64_t inline _mul128(int64_t a, int64_t b, int64_t* h) {
	uint64_t rhi;
	uint64_t rlo;
	__asm__("imulq  %[b];" :"=d"(rhi), "=a"(rlo) : "1"(a), [b]"rm"(b));
	*h = rhi;
	return rlo;
}

static uint64_t inline _udiv128(uint64_t hi, uint64_t lo, uint64_t d, uint64_t* r) {
	uint64_t q;
	uint64_t _r;
	__asm__("divq  %[d];" :"=d"(_r), "=a"(q) : "d"(hi), "a"(lo), [d]"rm"(d));
	*r = _r;
	return q;
}

static uint64_t inline __rdtsc() {
	uint32_t h;
	uint32_t l;
	__asm__("rdtsc;" :"=d"(h), "=a"(l));
	return (uint64_t)h << 32 | (uint64_t)l;
}

#define __shiftright128(a,b,n) ((a)>>(n))|((b)<<(64-(n)))
#define __shiftleft128(a,b,n) ((b)<<(n))|((a)>>(64-(n)))


#define _subborrow_u64(a,b,c,d) __builtin_ia32_sbb_u64(a,b,c,(long long unsigned int*)d);
#define _addcarry_u64(a,b,c,d) __builtin_ia32_addcarryx_u64(a,b,c,(long long unsigned int*)d);
#define _byteswap_uint64 __builtin_bswap64
#define LZC(x) __builtin_clzll(x)
#define TZC(x) __builtin_ctzll(x)

static void inline imm_mul(uint64_t* x, uint64_t y, uint64_t* dst, uint64_t* carryH) {
	unsigned char c = 0;
	uint64_t h, carry;
	dst[0] = _umul128(x[0], y, &h); carry = h;
	c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
	c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
	c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;
	c = _addcarry_u64(c, _umul128(x[4], y, &h), carry, dst + 4); carry = h;

	* carryH = carry;
}

static void inline imm_imul(uint64_t* x, uint64_t y, uint64_t* dst, uint64_t* carryH) {
	unsigned char c = 0;
	uint64_t h, carry;
	dst[0] = _umul128(x[0], y, &h); carry = h;
	c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
	c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
	c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;

	c = _addcarry_u64(c, _mul128(x[NB64BLOCK - 1], y, (int64_t*)&h), carry, dst + NB64BLOCK - 1); carry = h;
	*carryH = carry;
}

static void inline imm_umul(uint64_t* x, uint64_t y, uint64_t* dst) {
	// Assume that x[NB64BLOCK-1] is 0
	unsigned char c = 0;
	uint64_t h, carry;
	dst[0] = _umul128(x[0], y, &h); carry = h;
	c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
	c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
	c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;

	_addcarry_u64(c, 0ULL, carry, dst + (NB64BLOCK - 1));
}

static void inline shiftR(unsigned char n, uint64_t* d) {
	d[0] = __shiftright128(d[0], d[1], n);
	d[1] = __shiftright128(d[1], d[2], n);
	d[2] = __shiftright128(d[2], d[3], n);
	d[3] = __shiftright128(d[3], d[4], n);

	d[NB64BLOCK - 1] = ((int64_t)d[NB64BLOCK - 1]) >> n;
}

static void inline shiftR(unsigned char n, uint64_t* d, uint64_t h) {
	d[0] = __shiftright128(d[0], d[1], n);
	d[1] = __shiftright128(d[1], d[2], n);
	d[2] = __shiftright128(d[2], d[3], n);
	d[3] = __shiftright128(d[3], d[4], n);

	d[NB64BLOCK - 1] = __shiftright128(d[NB64BLOCK - 1], h, n);
}

//----------------- INT.CPP
Int _ONE((uint64_t)1);

Int::Int() {}

Int::Int(Int* a) {	if (a) Set(a);	else CLEAR();  }

Int::Int(int64_t i64) {
	if (i64 < 0) {	memset(bits64, 0xFF, NB64BLOCK * 8); }	else {	CLEAR(); }
	bits64[0] = i64;
}

Int::Int(uint64_t u64) {	CLEAR();	bits64[0] = u64; }
// ------------------------------------------------
void Int::CLEAR() { memset(bits64, 0, NB64BLOCK * 8); }
// ------------------------------------------------
void Int::Set(Int* a) {
	for (int i = 0; i < NB64BLOCK; i++){ 
		bits64[i] = a->bits64[i]; }	
}
// ------------------------------------------------ 
void Int::Add(Int* a) {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, bits64[4], a->bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Add(uint64_t a) {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a, bits64 + 0);
	c = _addcarry_u64(c, bits64[1], 0, bits64 + 1);
	c = _addcarry_u64(c, bits64[2], 0, bits64 + 2);
	c = _addcarry_u64(c, bits64[3], 0, bits64 + 3);
	c = _addcarry_u64(c, bits64[4], 0, bits64 + 4);
}
// ------------------------------------------------
void Int::AddOne() {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], 1, bits64 + 0);
	c = _addcarry_u64(c, bits64[1], 0, bits64 + 1);
	c = _addcarry_u64(c, bits64[2], 0, bits64 + 2);
	c = _addcarry_u64(c, bits64[3], 0, bits64 + 3);
	c = _addcarry_u64(c, bits64[4], 0, bits64 + 4);
}
// ------------------------------------------------
void Int::Add(Int* a, Int* b) {
	unsigned char c = 0;
	c = _addcarry_u64(c, b->bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, b->bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, b->bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, b->bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, b->bits64[4], a->bits64[4], bits64 + 4);
}
// ------------------------------------------------
uint64_t Int::AddCh(Int* a, uint64_t ca, Int* b, uint64_t cb) {
	uint64_t carry;
	unsigned char c = 0;
	c = _addcarry_u64(c, a->bits64[0], b->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, a->bits64[1], b->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, a->bits64[2], b->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, a->bits64[3], b->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, a->bits64[4], b->bits64[4], bits64 + 4);

	_addcarry_u64(c, ca, cb, &carry);
	return carry;
}

uint64_t Int::AddCh(Int* a, uint64_t ca) {
	uint64_t carry;
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, bits64[4], a->bits64[4], bits64 + 4);

	_addcarry_u64(c, ca, 0, &carry);
	return carry;
}
// ------------------------------------------------
uint64_t Int::AddC(Int* a) {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, bits64[4], a->bits64[4], bits64 + 4);

	return c;
}
// ------------------------------------------------
void Int::AddAndShift(Int* a, Int* b, uint64_t cH) {
	unsigned char c = 0;
	c = _addcarry_u64(c, b->bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, b->bits64[1], a->bits64[1], bits64 + 0);
	c = _addcarry_u64(c, b->bits64[2], a->bits64[2], bits64 + 1);
	c = _addcarry_u64(c, b->bits64[3], a->bits64[3], bits64 + 2);
	c = _addcarry_u64(c, b->bits64[4], a->bits64[4], bits64 + 3);

	bits64[NB64BLOCK - 1] = c + cH;
}
// ------------------------------------------------
void Int::MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22, uint64_t* cu, uint64_t* cv) {
	Int t1, t2, t3, t4;
	uint64_t c1, c2, c3, c4;
	c1 = t1.IMult(u, _11);
	c2 = t2.IMult(v, _12);
	c3 = t3.IMult(u, _21);
	c4 = t4.IMult(v, _22);
	*cu = u->AddCh(&t1, c1, &t2, c2);
	*cv = v->AddCh(&t3, c3, &t4, c4);
}

void Int::MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22) {
	Int t1, t2, t3, t4;
	t1.IMult(u, _11);
	t2.IMult(v, _12);
	t3.IMult(u, _21);
	t4.IMult(v, _22);
	u->Add(&t1, &t2);
	v->Add(&t3, &t4);
}

// ------------------------------------------------
bool Int::IsGreaterOrEqual(Int* a) {
	Int p;
	p.Sub(this, a);
	return p.IsPositive();
}
// ------------------------------------------------
bool Int::IsEqual(Int* a) {
	return
		(bits64[4] == a->bits64[4]) &&
		(bits64[3] == a->bits64[3]) &&
		(bits64[2] == a->bits64[2]) &&
		(bits64[1] == a->bits64[1]) &&
		(bits64[0] == a->bits64[0]);
}
// ------------------------------------------------
bool Int::IsOne() {	return IsEqual(&_ONE); }
// ------------------------------------------------
bool Int::IsZero() { return (bits64[4] | bits64[3] | bits64[2] | bits64[1] | bits64[0]) == 0; }
// ------------------------------------------------
void Int::SetInt32(uint32_t value) {	CLEAR();	bits[0] = value;	}
// ------------------------------------------------
unsigned char Int::GetByte(int n) {
	unsigned char* bbPtr = (unsigned char*)bits;
	return bbPtr[n];
}
// ------------------------------------------------
void Int::Get32Bytes(unsigned char* buff) {
	uint64_t* ptr = (uint64_t*)buff;
	ptr[3] = _byteswap_uint64(bits64[0]);
	ptr[2] = _byteswap_uint64(bits64[1]);
	ptr[1] = _byteswap_uint64(bits64[2]);
	ptr[0] = _byteswap_uint64(bits64[3]);
}
// ------------------------------------------------
void Int::Sub(Int* a) {
	unsigned char c = 0;
	c = _subborrow_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _subborrow_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _subborrow_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _subborrow_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _subborrow_u64(c, bits64[4], a->bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Sub(Int* a, Int* b) {
	unsigned char c = 0;
	c = _subborrow_u64(c, a->bits64[0], b->bits64[0], bits64 + 0);
	c = _subborrow_u64(c, a->bits64[1], b->bits64[1], bits64 + 1);
	c = _subborrow_u64(c, a->bits64[2], b->bits64[2], bits64 + 2);
	c = _subborrow_u64(c, a->bits64[3], b->bits64[3], bits64 + 3);
	c = _subborrow_u64(c, a->bits64[4], b->bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Sub(uint64_t a) {
	unsigned char c = 0;
	c = _subborrow_u64(c, bits64[0], a, bits64 + 0);
	c = _subborrow_u64(c, bits64[1], 0, bits64 + 1);
	c = _subborrow_u64(c, bits64[2], 0, bits64 + 2);
	c = _subborrow_u64(c, bits64[3], 0, bits64 + 3);
	c = _subborrow_u64(c, bits64[4], 0, bits64 + 4);
}
// ------------------------------------------------
bool Int::IsPositive() {	return (int64_t)(bits64[NB64BLOCK - 1]) >= 0;  }
// ------------------------------------------------
bool Int::IsNegative() { return (int64_t)(bits64[NB64BLOCK - 1]) < 0;   }
// ------------------------------------------------
bool Int::IsEven() {	return (bits[0] & 0x1) == 0;  }
// ------------------------------------------------
void Int::Neg() {
	unsigned char c = 0;
	c = _subborrow_u64(c, 0, bits64[0], bits64 + 0);
	c = _subborrow_u64(c, 0, bits64[1], bits64 + 1);
	c = _subborrow_u64(c, 0, bits64[2], bits64 + 2);
	c = _subborrow_u64(c, 0, bits64[3], bits64 + 3);
	c = _subborrow_u64(c, 0, bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Mult(Int* a) {	Int b(this);	Mult(a, &b);	}
// ------------------------------------------------
uint64_t Int::Mult(uint64_t a) {
	uint64_t carry;
	imm_mul(bits64, a, bits64, &carry);
	return carry;
}
// ------------------------------------------------
uint64_t Int::IMult(Int* a, int64_t b) {
	uint64_t carry;
	// Make b positive
	if (b < 0LL) {
		unsigned char c = 0;
		c = _subborrow_u64(c, 0, a->bits64[0], bits64 + 0);
		c = _subborrow_u64(c, 0, a->bits64[1], bits64 + 1);
		c = _subborrow_u64(c, 0, a->bits64[2], bits64 + 2);
		c = _subborrow_u64(c, 0, a->bits64[3], bits64 + 3);
		c = _subborrow_u64(c, 0, a->bits64[4], bits64 + 4);

		b = -b;
	} else {	Set(a);	}

	imm_imul(bits64, b, bits64, &carry);
	return carry;
}
// ------------------------------------------------
uint64_t Int::Mult(Int* a, uint64_t b) {
	uint64_t carry;
	imm_mul(a->bits64, b, bits64, &carry);
	return carry;
}
// ------------------------------------------------
void Int::Mult(Int* a, Int* b) {

	unsigned char c = 0;
	uint64_t h;
	uint64_t pr = 0;
	uint64_t carryh = 0;
	uint64_t carryl = 0;

	bits64[0] = _umul128(a->bits64[0], b->bits64[0], &pr);

	for (int i = 1; i < NB64BLOCK; i++) {
		for (int j = 0; j <= i; j++) {
			c = _addcarry_u64(c, _umul128(a->bits64[j], b->bits64[i - j], &h), pr, &pr);
			c = _addcarry_u64(c, carryl, h, &carryl);
			c = _addcarry_u64(c, carryh, 0, &carryh);
		}
		bits64[i] = pr;
		pr = carryl;
		carryl = carryh;
		carryh = 0;
	}
}
// ------------------------------------------------
int Int::GetSize() {
	int i = NB32BLOCK - 1;
	while (i > 0 && bits[i] == 0) i--;
	return i + 1;
}
// ------------------------------------------------
void Int::SetBase10(const char* value) {

	CLEAR();
	Int pw((uint64_t)1);
	Int c;
	int lgth = (int)strlen(value);
	for (int i = lgth - 1; i >= 0; i--) {
		uint32_t id = (uint32_t)(value[i] - '0');
		c.Set(&pw);
		c.Mult(id);
		Add(&c);
		pw.Mult(10);
	}
}

// ------------------------------------------------
void  Int::SetBase16(const char* value) {	SetBaseN(16, "0123456789ABCDEF", value); }
// ------------------------------------------------
std::string Int::GetBase10() {	return GetBaseN(10, "0123456789"); }
// ------------------------------------------------
std::string Int::GetBase16() {	return GetBaseN(16, "0123456789ABCDEF"); }
// ------------------------------------------------
void  Int::SetBaseN(int n, const char* charset, const char* value) {

	CLEAR();

	Int pw((uint64_t)1);
	Int nb((uint64_t)n);
	Int c;

	int lgth = (int)strlen(value);
	for (int i = lgth - 1; i >= 0; i--) {
		const char* p = strchr(charset, toupper(value[i]));
		if (!p) { printf("Invalid charset !!\n");	return;	}
		int id = (int)(p - charset);
		c.SetInt32(id);
		c.Mult(&pw);
		Add(&c);
		pw.Mult(&nb);

	}
}
// ------------------------------------------------
std::string Int::GetBaseN(int n, char* charset) {

	std::string ret;

	Int N(this);
	int isNegative = N.IsNegative();
	if (isNegative) N.Neg();

	// TODO: compute max digit
	unsigned char digits[1024];
	memset(digits, 0, sizeof(digits));

	int digitslen = 1;
	for (int i = 0; i < NB64BLOCK * 8; i++) {
		unsigned int carry = N.GetByte(NB64BLOCK * 8 - i - 1);
		for (int j = 0; j < digitslen; j++) {
			carry += (unsigned int)(digits[j]) << 8;
			digits[j] = (unsigned char)(carry % n);
			carry /= n;
		}
		while (carry > 0) {
			digits[digitslen++] = (unsigned char)(carry % n);
			carry /= n;
		}
	}

	// reverse
	if (isNegative)
		ret.push_back('-');

	for (int i = 0; i < digitslen; i++)
		ret.push_back(charset[digits[digitslen - 1 - i]]);

	if (ret.length() == 0)
		ret.push_back('0');

	return ret;
}
//============== INT.CPP - end ================================================================================

//============== INT->INTMOD.CPP - start ================================================================================
#include <emmintrin.h>

static Int     _P;       // Field characteristic
static Int     _R;       // Montgomery multiplication R
static Int     _R2;      // Montgomery multiplication R2
static Int     _R3;      // Montgomery multiplication R3
static Int     _R4;      // Montgomery multiplication R4
static int32_t  Msize;    // Montgomery mult size
static uint64_t MM64;     // 64bits lsb negative inverse of P
#define MSK62  0x3FFFFFFFFFFFFFFF

extern Int _ONE;

// ------------------------------------------------
void Int::ModAdd(Int* a) {
	Int p;

	Add(a);
	p.Sub(this, &_P);
	if (p.IsPositive())
		Set(&p);
}
// ------------------------------------------------
void Int::ModAdd(Int* a, Int* b) {
	Int p;
	Add(a, b);
	p.Sub(this, &_P);
	if (p.IsPositive()){Set(&p);}
		
}
// ------------------------------------------------
void Int::ModSub(Int* a) {	Sub(a); 	if (IsNegative()){ Add(&_P); }	}
// ------------------------------------------------
void Int::ModSub(Int* a, Int* b) {	Sub(a, b);	if (IsNegative()){ Add(&_P); }   }
// ------------------------------------------------
void Int::ModNeg() {	Neg();	Add(&_P);	}
// ------------------------------------------------
void Int::DivStep62(Int* u, Int* v, int64_t* eta, int* pos, int64_t* uu, int64_t* uv, int64_t* vu, int64_t* vv) {
	int bitCount;
	uint64_t u0 = u->bits64[0];
	uint64_t v0 = v->bits64[0];

#define SWAP(tmp,x,y) tmp = x; x = y; y = tmp;
	uint64_t uh, vh, w, x;
	unsigned char c = 0;

	// Extract 64 MSB of u and v
	// u and v must be positive

	while (*pos >= 1 && (u->bits64[*pos] | v->bits64[*pos]) == 0) (*pos)--;
	if (*pos == 0) {	uh = u->bits64[0];	vh = v->bits64[0];	}
	else {
		uint64_t s = LZC(u->bits64[*pos] | v->bits64[*pos]);
		if (s == 0) {	uh = u->bits64[*pos];	vh = v->bits64[*pos];	}
		else {
			uh = __shiftleft128(u->bits64[*pos - 1], u->bits64[*pos], (uint8_t)s);
			vh = __shiftleft128(v->bits64[*pos - 1], v->bits64[*pos], (uint8_t)s);
		}
	}

	bitCount = 62;
	__m128i _u, _v, _t;

	((int64_t*)&_u)[0] = 1;
	((int64_t*)&_u)[1] = 0;
	((int64_t*)&_v)[0] = 0;
	((int64_t*)&_v)[1] = 1;

	while (true) {
		// Use a sentinel bit to count zeros only up to bitCount
		uint64_t zeros = TZC(v0 | 1ULL << bitCount);
		vh >>= zeros;
		v0 >>= zeros;
		_u = _mm_slli_epi64(_u, (int)zeros);
		bitCount -= (int)zeros;

		if (bitCount <= 0) {	break;	}

		if (vh < uh) {
			SWAP(w, uh, vh);
			SWAP(x, u0, v0);
			SWAP(_t, _u, _v);
		}

		vh -= uh;
		v0 -= u0;
		_v = _mm_sub_epi64(_v, _u);
	}

	*uu = ((int64_t*)&_u)[0];
	*uv = ((int64_t*)&_u)[1];
	*vu = ((int64_t*)&_v)[0];
	*vv = ((int64_t*)&_v)[1];
}

// ------------------------------------------------
void Int::ModInv() {

	Int u(&_P);
	Int v(this);
	Int r((int64_t)0);
	Int s((int64_t)1);

	// Delayed right shift 62bits
	Int r0_P;
	Int s0_P;

	int64_t  eta = -1;
	int64_t uu, uv, vu, vv;
	uint64_t carryS, carryR;
	int pos = NB64BLOCK - 1;
	while (pos >= 1 && (u.bits64[pos] | v.bits64[pos]) == 0) pos--;

	while (!v.IsZero()) { 
		DivStep62(&u, &v, &eta, &pos, &uu, &uv, &vu, &vv);
		// Now update BigInt variables
		MatrixVecMul(&u, &v, uu, uv, vu, vv);

		// Make u,v positive
		// Required only for Pornin's method
		if (u.IsNegative()) {
			u.Neg();
			uu = -uu;
			uv = -uv;
		}
		if (v.IsNegative()) {
			v.Neg();
			vu = -vu;
			vv = -vv;
		}

		MatrixVecMul(&r, &s, uu, uv, vu, vv, &carryR, &carryS);

		// Compute multiple of P to add to s and r to make them multiple of 2^62
		uint64_t r0 = (r.bits64[0] * MM64) & MSK62;
		uint64_t s0 = (s.bits64[0] * MM64) & MSK62;
		r0_P.Mult(&_P, r0);
		s0_P.Mult(&_P, s0);
		carryR = r.AddCh(&r0_P, carryR);
		carryS = s.AddCh(&s0_P, carryS);

		// Right shift all variables by 62bits
		shiftR(62, u.bits64);
		shiftR(62, v.bits64);
		shiftR(62, r.bits64, carryR);
		shiftR(62, s.bits64, carryS);

	}

	if (!u.IsOne()) { CLEAR(); return;	}		// No inverse
	while (r.IsNegative()){ r.Add(&_P); }
	while (r.IsGreaterOrEqual(&_P)){ r.Sub(&_P); }
	Set(&r);
}
// ------------------------------------------------
void Int::ModMul(Int* a, Int* b) {
	Int p;
	p.MontgomeryMult(a, b);
	MontgomeryMult(&_R2, &p);
}
// ------------------------------------------------
void Int::SetupField(Int* n, Int* R, Int* R2, Int* R3, Int* R4) {

	// Size in number of 32bit word
	int nSize = n->GetSize();

	// Last digit inversions (Newton's iteration)
	{
		int64_t x, t;
		x = t = (int64_t)n->bits64[0];
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		MM64 = (uint64_t)(-x);
	}
	_P.Set(n);

	// Size of Montgomery mult (64bits digit)
	Msize = nSize / 2;

	// Compute few power of R
	// R = 2^(64*Msize) mod n
	Int Ri;
	Ri.MontgomeryMult(&_ONE, &_ONE); // Ri = R^-1
	_R.Set(&Ri);                     // R  = R^-1
	_R2.MontgomeryMult(&Ri, &_ONE);  // R2 = R^-2
	_R3.MontgomeryMult(&Ri, &Ri);    // R3 = R^-3
	_R4.MontgomeryMult(&_R3, &_ONE); // R4 = R^-4

	_R.ModInv();                     // R  = R
	_R2.ModInv();                    // R2 = R^2
	_R3.ModInv();                    // R3 = R^3
	_R4.ModInv();                    // R4 = R^4

	if (R){ R->Set(&_R); }
	if (R2){R2->Set(&_R2);}
	if (R3){R3->Set(&_R3);}
	if (R4){R4->Set(&_R4);}
}
// ------------------------------------------------
void Int::MontgomeryMult(Int* a, Int* b) {
	Int pr;
	Int p;
	uint64_t ML;
	uint64_t c;

	// i = 0
	imm_umul(a->bits64, b->bits64[0], pr.bits64);
	ML = pr.bits64[0] * MM64;
	imm_umul(_P.bits64, ML, p.bits64);
	c = pr.AddC(&p);
	memcpy(bits64, pr.bits64 + 1, 8 * (NB64BLOCK - 1));
	bits64[NB64BLOCK - 1] = c;

	for (int i = 1; i < Msize; i++) {
		imm_umul(a->bits64, b->bits64[i], pr.bits64);
		ML = (pr.bits64[0] + bits64[0]) * MM64;
		imm_umul(_P.bits64, ML, p.bits64);
		c = pr.AddC(&p);
		AddAndShift(this, &pr, c);
	}

	p.Sub(this, &_P);
	if (p.IsPositive()){ Set(&p); }
}

// SecpK1 specific section -----------------------------------------------------------------------------
void Int::ModMulK1(Int* a, Int* b) 
{
	unsigned char c;
	uint64_t ah, al;
	uint64_t t[NB64BLOCK];
	uint64_t r512[8];
	r512[5] = 0;
	r512[6] = 0;
	r512[7] = 0;

	// 256*256 multiplier
	imm_umul(a->bits64, b->bits64[0], r512);
	imm_umul(a->bits64, b->bits64[1], t);
	c = _addcarry_u64(0, r512[1], t[0], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[1], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[2], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[3], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[4], r512 + 5);
	imm_umul(a->bits64, b->bits64[2], t);
	c = _addcarry_u64(0, r512[2], t[0], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[1], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[2], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[3], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[4], r512 + 6);
	imm_umul(a->bits64, b->bits64[3], t);
	c = _addcarry_u64(0, r512[3], t[0], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[1], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[2], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[3], r512 + 6);
	c = _addcarry_u64(c, r512[7], t[4], r512 + 7);

	// Reduce from 512 to 320 
	imm_umul(r512 + 4, 0x1000003D1ULL, t);
	c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
	c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

	// Reduce from 320 to 256 
	// No overflow possible here t[4]+c<=0x1000003D1ULL
	al = _umul128(t[4] + c, 0x1000003D1ULL, &ah);
	c = _addcarry_u64(0, r512[0], al, bits64 + 0);
	c = _addcarry_u64(c, r512[1], ah, bits64 + 1);
	c = _addcarry_u64(c, r512[2], 0ULL, bits64 + 2);
	c = _addcarry_u64(c, r512[3], 0ULL, bits64 + 3);

	// Probability of carry here or that this>P is very very unlikely
	bits64[4] = 0;
}
// ------------------------------------------------
void Int::ModMulK1(Int* a) 
{
	unsigned char c;
	uint64_t ah, al;
	uint64_t t[NB64BLOCK];
	uint64_t r512[8];
	r512[5] = 0;
	r512[6] = 0;
	r512[7] = 0;

	// 256*256 multiplier
	imm_umul(a->bits64, bits64[0], r512);
	imm_umul(a->bits64, bits64[1], t);
	c = _addcarry_u64(0, r512[1], t[0], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[1], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[2], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[3], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[4], r512 + 5);
	imm_umul(a->bits64, bits64[2], t);
	c = _addcarry_u64(0, r512[2], t[0], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[1], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[2], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[3], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[4], r512 + 6);
	imm_umul(a->bits64, bits64[3], t);
	c = _addcarry_u64(0, r512[3], t[0], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[1], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[2], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[3], r512 + 6);
	c = _addcarry_u64(c, r512[7], t[4], r512 + 7);

	// Reduce from 512 to 320 
	imm_umul(r512 + 4, 0x1000003D1ULL, t);
	c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
	c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

	// Reduce from 320 to 256 
	// No overflow possible here t[4]+c<=0x1000003D1ULL
	al = _umul128(t[4] + c, 0x1000003D1ULL, &ah);
	c = _addcarry_u64(0, r512[0], al, bits64 + 0);
	c = _addcarry_u64(c, r512[1], ah, bits64 + 1);
	c = _addcarry_u64(c, r512[2], 0, bits64 + 2);
	c = _addcarry_u64(c, r512[3], 0, bits64 + 3);
	// Probability of carry here or that this>P is very very unlikely
	bits64[4] = 0;
}
// ------------------------------------------------
void Int::ModSquareK1(Int* a) 
{
	unsigned char c;
	uint64_t t[NB64BLOCK], SL, SH, r512[8], t1, t2;

	//k=0
	r512[0] = _umul128(a->bits64[0], a->bits64[0], &t[1]); 

	//k=1
	t[3] = _umul128(a->bits64[0], a->bits64[1], &t[4]);
	c = _addcarry_u64(0, t[3], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], t[4], &t[4]);
	c = _addcarry_u64(c, 0, 0, &t1);
	c = _addcarry_u64(0, t[1], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], 0, &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	r512[1] = t[3];

	//k=2
	t[0] = _umul128(a->bits64[0], a->bits64[2], &t[1]);
	c = _addcarry_u64(0, t[0], t[0], &t[0]);
	c = _addcarry_u64(c, t[1], t[1], &t[1]);
	c = _addcarry_u64(c, 0, 0, &t2);

	SL = _umul128(a->bits64[1], a->bits64[1], &SH);
	c = _addcarry_u64(0, t[0], SL, &t[0]);
	c = _addcarry_u64(c, t[1], SH, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	c = _addcarry_u64(0, t[0], t[4], &t[0]);
	c = _addcarry_u64(c, t[1], t1, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	r512[2] = t[0];

	//k=3
	t[3] = _umul128(a->bits64[0], a->bits64[3], &t[4]);
	SL = _umul128(a->bits64[1], a->bits64[2], &SH);

	c = _addcarry_u64(0, t[3], SL, &t[3]);
	c = _addcarry_u64(c, t[4], SH, &t[4]);
	c = _addcarry_u64(c, 0, 0, &t1);
	t1 += t1;
	c = _addcarry_u64(0, t[3], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], t[4], &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	c = _addcarry_u64(0, t[3], t[1], &t[3]);
	c = _addcarry_u64(c, t[4], t2, &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	r512[3] = t[3];

	//k=4
	t[0] = _umul128(a->bits64[1], a->bits64[3], &t[1]);
	c = _addcarry_u64(0, t[0], t[0], &t[0]);
	c = _addcarry_u64(c, t[1], t[1], &t[1]);
	c = _addcarry_u64(c, 0, 0, &t2);

	SL = _umul128(a->bits64[2], a->bits64[2], &SH);
	c = _addcarry_u64(0, t[0], SL, &t[0]);
	c = _addcarry_u64(c, t[1], SH, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	c = _addcarry_u64(0, t[0], t[4], &t[0]);
	c = _addcarry_u64(c, t[1], t1, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	r512[4] = t[0];

	//k=5
	t[3] = _umul128(a->bits64[2], a->bits64[3], &t[4]);
	c = _addcarry_u64(0, t[3], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], t[4], &t[4]);
	c = _addcarry_u64(c, 0, 0, &t1);
	c = _addcarry_u64(0, t[3], t[1], &t[3]);
	c = _addcarry_u64(c, t[4], t2, &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	r512[5] = t[3];

	//k=6
	t[0] = _umul128(a->bits64[3], a->bits64[3], &t[1]);
	c = _addcarry_u64(0, t[0], t[4], &t[0]);
	c = _addcarry_u64(c, t[1], t1, &t[1]);
	r512[6] = t[0];

	//k=7
	r512[7] = t[1];

	// Reduce from 512 to 320 
	imm_umul(r512 + 4, 0x1000003D1ULL, t);
	c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
	c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

	// Reduce from 320 to 256 
	// No overflow possible here t[4]+c<=0x1000003D1ULL
	SL = _umul128(t[4] + c, 0x1000003D1ULL, &SH);
	c = _addcarry_u64(0, r512[0], SL, bits64 + 0);
	c = _addcarry_u64(c, r512[1], SH, bits64 + 1);
	c = _addcarry_u64(c, r512[2], 0, bits64 + 2);
	c = _addcarry_u64(c, r512[3], 0, bits64 + 3);
	// Probability of carry here or that this>P is very very unlikely
	bits64[4] = 0;
}

//============== INT->INTMOD.CPP - end ================================================================================
//============== POINT - start ================================================================================
class Point
{
	public:
		Point();
		Point(const Point &p);
		~Point();

		void Clear();
		void Reduce();
		std::string toString();

		Int x, y, z;
};

//------------------
Point::Point(){}
Point::~Point(){}
//------------------
Point::Point(const Point &p){
    x.Set((Int *)&p.x);
    y.Set((Int *)&p.y);
    z.Set((Int *)&p.z);
}
//------------------
void Point::Clear(){
    x.SetInt32(0);
    y.SetInt32(0);
    z.SetInt32(0);
}
//------------------
void Point::Reduce(){
    Int i(&z);
    i.ModInv();
    x.ModMul(&x, &i);
    y.ModMul(&y, &i);
    z.SetInt32(1);
}
//------------------
std::string Point::toString()
{
    std::string ret;
    ret  = "\nX=" + x.GetBase16() + "\n";
    ret += "Y=" + y.GetBase16() + "\n";
    ret += "Z=" + z.GetBase16() + "\n";
    return ret;
}
//============== POINT - end ================================================================================

//============== CSHA256 - start ================================================================================

class CSHA256
	{
	private:
		uint32_t s[8];
		unsigned char buf[64];
		uint64_t bytes;

	public:
		static const size_t OUTPUT_SIZE = 32;

		CSHA256();
		void Write(const unsigned char* data, size_t len);
		void Finalize(unsigned char hash[OUTPUT_SIZE]);

	};

CSHA256::CSHA256() {
	bytes = 0;
	s[0] = 0x6a09e667ul;
	s[1] = 0xbb67ae85ul;
	s[2] = 0x3c6ef372ul;
	s[3] = 0xa54ff53aul;
	s[4] = 0x510e527ful;
	s[5] = 0x9b05688cul;
	s[6] = 0x1f83d9abul;
	s[7] = 0x5be0cd19ul;
}

void CSHA256::Write(const unsigned char* data, size_t len)
{
const unsigned char* end = data + len;
size_t bufsize = bytes % 64;
if (bufsize && bufsize + len >= 64) {
	// Fill the buffer, and process it.
	memcpy(buf + bufsize, data, 64 - bufsize);
	bytes += 64 - bufsize;
	data += 64 - bufsize;
	_sha256::Transform(s, buf);
	bufsize = 0;
}
while (end >= data + 64) {
	// Process full chunks directly from the source.
	_sha256::Transform(s, data);
	bytes += 64;
	data += 64;
}
if (end > data) {
	// Fill the buffer with what remains.
	memcpy(buf + bufsize, data, end - data);
	bytes += end - data;
}
}

void CSHA256::Finalize(unsigned char hash[OUTPUT_SIZE])
{
	unsigned char sizedesc[8];
	WRITEBE64(sizedesc, bytes << 3);
	Write(_sha256::pad, 1 + ((119 - (bytes % 64)) % 64));
	Write(sizedesc, 8);
	WRITEBE32(hash, s[0]);
	WRITEBE32(hash + 4, s[1]);
	WRITEBE32(hash + 8, s[2]);
	WRITEBE32(hash + 12, s[3]);
	WRITEBE32(hash + 16, s[4]);
	WRITEBE32(hash + 20, s[5]);
	WRITEBE32(hash + 24, s[6]);
	WRITEBE32(hash + 28, s[7]);
}

void sha256(unsigned char *input, int length, unsigned char *digest) {

	CSHA256 sha;
	sha.Write(input, length);
	sha.Finalize(digest);
}
//============== CSHA256 - end ================================================================================
//============== SECP256K - start ================================================================================
class Secp256K1
{
	public:
		Secp256K1();
		~Secp256K1();
		void Init();
		Point ComputePublicKey(Int* privKey);
		void GetHash160(int type, bool isCompressed, Point& pubKey, unsigned char* hash);

		Point Add2(Point& p1, Point& p2);
		Point AddDirect(Point& p1, Point& p2);
		Point DoubleDirect(Point& p);
		Point G;                 // Generator

	private:
		uint8_t GetByte(std::string& str, int idx);
		Point GTable[256 * 32];     // Generator table
};
//-------------------------------------------------------
Secp256K1::Secp256K1(){}
Secp256K1::~Secp256K1(){}

//-------------------------------------------------------
void Secp256K1::Init() 
{
	// Prime for the finite field
	Int P;
	P.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F");

	// Set up field
	Int::SetupField(&P);

	// Generator point // Điểm cơ sở = 04 79BE667E F9DCBBAC 55A06295 CE870B07 029BFCDB 2DCE28D9 59F2815B 16F81798 483ADA77 26A3C465 5DA4FBFC 0E1108A8 FD17B448 A6855419 9C47D08F FB10D4B8
	G.x.SetBase16("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798");
	G.y.SetBase16("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8");
	G.z.SetInt32(1);

	// Compute Generator table
	Point N(G);
	for (int i = 0; i < 32; i++) {
		GTable[i * 256] = N;
		N = DoubleDirect(N);
		for (int j = 1; j < 255; j++) {
			GTable[i * 256 + j] = N;
			N = AddDirect(N, GTable[i * 256]);
		}
		GTable[i * 256 + 255] = N; // Dummy point for check function
	}
}
//-------------------------------------------------------

Point Secp256K1::ComputePublicKey(Int* privKey)
{
	int i = 0;
	uint8_t b;
	Point Q;

	Q.Clear(); //Q ( x=0 , y=0 , z=0 )

	// Search first significant byte
	for (i = 0; i < 32; i++) {	
		b = privKey->GetByte(i);
		if (b){ break; }	}
		
	Q = GTable[256 * i + (b - 1)];
	i++;

	for (; i < 32; i++) {	
		b = privKey->GetByte(i);
		if (b){ Q = Add2(Q, GTable[256 * i + (b - 1)]); }	}

	Q.Reduce();

	return Q;
}
//-------------------------------------------------------

uint8_t Secp256K1::GetByte(std::string& str, int idx)
{
	char tmp[3];
	int  val;

	tmp[0] = str.data()[2 * idx];
	tmp[1] = str.data()[2 * idx + 1];
	tmp[2] = 0;

	if (sscanf(tmp, "%X", &val) != 1) {	printf("ParsePublicKeyHex: Error invalid public key specified (unexpected hexadecimal digit)\n");	exit(-1);	}

	return (uint8_t)val;
}
//-------------------------------------------------------

void Secp256K1::GetHash160(int type, bool isCompressed, Point& pubKey, unsigned char* hash)
{
	unsigned char shapk[64];

	switch (type) {
		case P2PKH:
		case BECH32:
			{
			unsigned char publicKeyBytes[128];

			if (!isCompressed) {
				// Full public key
				publicKeyBytes[0] = 0x4;
				pubKey.x.Get32Bytes(publicKeyBytes + 1);
				pubKey.y.Get32Bytes(publicKeyBytes + 33);
				sha256_65(publicKeyBytes, shapk);

			} else {
				// Compressed public key
				publicKeyBytes[0] = pubKey.y.IsEven() ? 0x2 : 0x3;
				pubKey.x.Get32Bytes(publicKeyBytes + 1);
				sha256_33(publicKeyBytes, shapk);
			}

			ripemd160_32(shapk, hash);
			}
		break;

		case P2SH:
			{
			// Redeem Script (1 to 1 P2SH)
			unsigned char script[64];

			script[0] = 0x00;  // OP_0
			script[1] = 0x14;  // PUSH 20 bytes
			GetHash160(P2PKH, isCompressed, pubKey, script + 2);

			sha256(script, 22, shapk);
			ripemd160_32(shapk, hash);
			}
		break;
	}

}

//-------------------------------------------------------

Point Secp256K1::AddDirect(Point& p1, Point& p2)
{
	Int _s, _p, dy, dx;
	Point r;
	r.z.SetInt32(1);

	dy.ModSub(&p2.y, &p1.y);
	dx.ModSub(&p2.x, &p1.x);
	dx.ModInv();
	_s.ModMulK1(&dy, &dx);    // s = (p2.y-p1.y)*inverse(p2.x-p1.x);

	_p.ModSquareK1(&_s);       // _p = pow2(s)

	r.x.ModSub(&_p, &p1.x);
	r.x.ModSub(&p2.x);       // rx = pow2(s) - p1.x - p2.x;

	r.y.ModSub(&p2.x, &r.x);
	r.y.ModMulK1(&_s);
	r.y.ModSub(&p2.y);       // ry = - p2.y - s*(ret.x-p2.x);

	return r;
}
//-------------------------------------------------------

Point Secp256K1::Add2(Point& p1, Point& p2)
{
	// P2.z = 1
	Int u, v, u1, v1, vs2, vs3, us2, a, us2w, vs2v2, vs3u2, _2vs2v2;
	Point r;

	u1.ModMulK1(&p2.y, &p1.z);
	v1.ModMulK1(&p2.x, &p1.z);
	u.ModSub(&u1, &p1.y);
	v.ModSub(&v1, &p1.x);
	us2.ModSquareK1(&u);
	vs2.ModSquareK1(&v);
	vs3.ModMulK1(&vs2, &v);
	us2w.ModMulK1(&us2, &p1.z);
	vs2v2.ModMulK1(&vs2, &p1.x);
	_2vs2v2.ModAdd(&vs2v2, &vs2v2);
	a.ModSub(&us2w, &vs3);
	a.ModSub(&_2vs2v2);

	r.x.ModMulK1(&v, &a);

	vs3u2.ModMulK1(&vs3, &p1.y);
	r.y.ModSub(&vs2v2, &a);
	r.y.ModMulK1(&r.y, &u);
	r.y.ModSub(&vs3u2);

	r.z.ModMulK1(&vs3, &p1.z);

	return r;
}
//-------------------------------------------------------
 
Point Secp256K1::DoubleDirect(Point& p)
{
	Int _s, _p, a;
	Point r;
	r.z.SetInt32(1);

	_s.ModMulK1(&p.x, &p.x);
	_p.ModAdd(&_s, &_s);
	_p.ModAdd(&_s);

	a.ModAdd(&p.y, &p.y);
	a.ModInv();
	_s.ModMulK1(&_p, &a);    // s = (3*pow2(p.x))*inverse(2*p.y);

	_p.ModMulK1(&_s, &_s);
	a.ModAdd(&p.x, &p.x);
	a.ModNeg();
	r.x.ModAdd(&a, &_p);   // rx = pow2(s) + neg(2*p.x);

	a.ModSub(&r.x, &p.x);

	_p.ModMulK1(&a, &_s);
	r.y.ModAdd(&_p, &p.y);
	r.y.ModNeg();           // ry = neg(p.y + s*(ret.x+neg(p.x)));

	return r;
}
//============== SECP256K - end ================================================================================
//============== BECH_32 - start ================================================================================
	static const char* charset = "qpzry9x8gf2tvdw0s3jn54khce6mua7l";
	
//----------------------------------------------------------
	static int convert_bits(uint8_t* out, size_t* outlen, int outbits, const uint8_t* in, size_t inlen, int inbits, int pad) {
	uint32_t val = 0;
	int bits = 0;
	uint32_t maxv = (((uint32_t)1) << outbits) - 1;
	while (inlen--) {
		val = (val << inbits) | *(in++);
		bits += inbits;
		while (bits >= outbits) {
		bits -= outbits;
		out[(*outlen)++] = (val >> bits) & maxv;
		}
	}
	if (pad) {
		if (bits) {
		out[(*outlen)++] = (val << (outbits - bits)) & maxv;
		}
	} else if (((val << (outbits - bits)) & maxv) || bits >= inbits) {
		return 0;
	}
	return 1;
	}

//----------------------------------------------------------
	uint32_t bech32_polymod_step(uint32_t pre) {
	uint8_t b = pre >> 25;
	return ((pre & 0x1FFFFFF) << 5) ^
		(-((b >> 0) & 1) & 0x3b6a57b2UL) ^
		(-((b >> 1) & 1) & 0x26508e6dUL) ^
		(-((b >> 2) & 1) & 0x1ea119faUL) ^
		(-((b >> 3) & 1) & 0x3d4233ddUL) ^
		(-((b >> 4) & 1) & 0x2a1462b3UL);
	}

//----------------------------------------------------------
	int bech32_encode(char *output, const char *hrp, const uint8_t *data, size_t data_len) {
	uint32_t chk = 1;
	size_t i = 0;
	while (hrp[i] != 0) {
		int ch = hrp[i];
		if (ch < 33 || ch > 126) {
		return 0;
		}

		if (ch >= 'A' && ch <= 'Z') return 0;
		chk = bech32_polymod_step(chk) ^ (ch >> 5);
		++i;
	}
	if (i + 7 + data_len > 90) return 0;
	chk = bech32_polymod_step(chk);
	while (*hrp != 0) {
		chk = bech32_polymod_step(chk) ^ (*hrp & 0x1f);
		*(output++) = *(hrp++);
	}
	*(output++) = '1';
	for (i = 0; i < data_len; ++i) {
		if (*data >> 5) return 0;
		chk = bech32_polymod_step(chk) ^ (*data);
		*(output++) = charset[*(data++)];
	}
	for (i = 0; i < 6; ++i) {
		chk = bech32_polymod_step(chk);
	}
	chk ^= 1;
	for (i = 0; i < 6; ++i) {
		*(output++) = charset[(chk >> ((5 - i) * 5)) & 0x1f];
	}
	*output = 0;
	return 1;
	}

//----------------------------------------------------------
int segwit_addr_encode(char *output, const char *hrp, int witver, const uint8_t *witprog, size_t witprog_len) {
	uint8_t data[65];
	size_t datalen = 0;
	if (witver > 16) return 0;
	if (witver == 0 && witprog_len != 20 && witprog_len != 32) return 0;
	if (witprog_len < 2 || witprog_len > 40) return 0;
	data[0] = witver;
	convert_bits(data + 1, &datalen, 5, witprog, witprog_len, 8, 1);
	++datalen;
	return bech32_encode(output, hrp, data, datalen);
}

//============== BECH_32 - end ================================================================================


//============== HIIU =============================================================================
//============== HIIU =============================================================================
//============== HIIU =============================================================================
//============== HIIU =============================================================================
//============== HIIU =============================================================================
//============== HIIU =============================================================================
//============== HIIU =============================================================================

class Bitcoin 
{
	public:
		Point privToPubkey(char* p_hex); //return POINT
		void privToHash160(int type, char* p_hex, uint32_t* _hash160, bool isCompressed);
		std::string privToAddr(int type, char* p_hex, bool isCompressed); // return STRING

		void pubkeyToHash160(int type, Point& pubKey, uint32_t* _hash160, bool isCompressed); 
		std::string pubkeyToAddr(int type, Point& pubKey, bool isCompressed); //return STRING 

		std::string hash160ToAddr(uint32_t* _hash160, bool isCompressed); //return STRING 
};  
//-------------------------------------------------------------------- 

Point Bitcoin::privToPubkey(char* p_hex)
{
	char* priv_hex = p_hex;

	Int priv; 
	priv.SetBase16(priv_hex); 
	 
	Secp256K1*	hiiu_secp = new Secp256K1();   
	hiiu_secp->Init();	

	Point pubKey = hiiu_secp->ComputePublicKey(&priv); //

	delete hiiu_secp;

	return pubKey;
}

//--------------------------------------------------------------------
 void Bitcoin::privToHash160(int type, char* p_hex, uint32_t* _hash160, bool isCompressed)
 {

	char* priv_hex = p_hex;   // char* priv_hex = "2832ed74f2b5e35ee"; 

	Int priv; 
	priv.SetBase16(priv_hex); 
	
	Secp256K1*	hiiu_secp = new Secp256K1();   
	hiiu_secp->Init();	

	Point pubKey = hiiu_secp->ComputePublicKey(&priv); //

	unsigned char address[25];

	switch (type) { 
		case P2PKH:
			break;	
		case P2SH:
			if (!isCompressed) {  std::cout<< "\n P2SH: Only compressed key "; exit(-1); }
			break;	
		case BECH32:
			if (!isCompressed) {std::cout<< "\n BECH32: Only compressed key "; exit(-1);	}
			break;
	}

	hiiu_secp->GetHash160(type, isCompressed, pubKey, address + 1);
	
	//print 
	printf("\n hash160 : ");	
	for (int i = 1; i < 21; i++){	
		printf("%.2x", address[i]);		
	}

	// 20 bytes --> uint32_t _hash160[5]
	memcpy(_hash160, address + 1, 4);
	memcpy(_hash160 + 1, address + 5, 4);
	memcpy(_hash160 + 2, address + 9, 4);
	memcpy(_hash160 + 3, address + 13, 4);
	memcpy(_hash160 + 4, address + 17, 4);

	delete hiiu_secp;
 }

//--------------------------------------------------------------------

std::string Bitcoin::privToAddr(int type, char* p_hex, bool isCompressed)
{
	char* priv_hex = p_hex;   // char* priv_hex = "2832ed74f2b5e35ee"; 

	Int priv; 
	priv.SetBase16(priv_hex); 
	
	Secp256K1*	hiiu_secp = new Secp256K1();   
	hiiu_secp->Init();	

	Point pubKey = hiiu_secp->ComputePublicKey(&priv); //

	unsigned char address[25];

	switch (type) { 
		case P2PKH:
			address[0] = 0x00;
			break;
		case P2SH:
			if (!isCompressed) {  return " P2SH: Only compressed key ";  }
			address[0] = 0x05;
			break;	
		case BECH32:
			if (!isCompressed) {	return " BECH32: Only compressed key ";	}
			hiiu_secp->GetHash160(type, isCompressed, pubKey, address + 1);

			char addr_bech32[128];
			segwit_addr_encode(addr_bech32, "bc", 0, address + 1, 20);
			return std::string(addr_bech32); 
	}

	hiiu_secp->GetHash160(type, isCompressed, pubKey, address + 1);	

	sha256_checksum(address, 21, address + 21);

	std::string addr = EncodeBase58(address, address + 25);

	delete hiiu_secp;
	return addr;
 }
//--------------------------------------------------------------------
void Bitcoin::pubkeyToHash160(int type, Point& pubKey, uint32_t* _hash160,  bool isCompressed)
{
	Secp256K1* hiiu_secp = new Secp256K1();   
	hiiu_secp->Init();	
	
	unsigned char address[25];

	switch (type) { 
		case P2PKH:
			break;	
		case P2SH:
			if (!isCompressed) {  std::cout<< "\n P2SH: Only compressed key "; exit(-1); }
			break;	
		case BECH32:
			if (!isCompressed) {std::cout<< "\n BECH32: Only compressed key "; exit(-1);	}
			break;
	}

	hiiu_secp->GetHash160(type, isCompressed, pubKey, address + 1);
	
	//print 
	printf("\n hash160 : ");
	for (int i = 1; i < 21; i++){	
		printf("%.2x", address[i]);			
	}

	// tách 20 bytes --> uint32_t _hash160[5]
	memcpy(_hash160, address + 1, 4);
	memcpy(_hash160 + 1, address + 5, 4);
	memcpy(_hash160 + 2, address + 9, 4);
	memcpy(_hash160 + 3, address + 13, 4);
	memcpy(_hash160 + 4, address + 17, 4);

	delete hiiu_secp;
} 
//--------------------------------------------------------------------
std::string Bitcoin::pubkeyToAddr(int type, Point& pubKey, bool isCompressed)
{
	Secp256K1*	hiiu_secp = new Secp256K1();   
	hiiu_secp->Init();	
	
	unsigned char address[25];

	switch (type) { 
		case P2PKH:
			address[0] = 0x00;
			break;			
		case P2SH:
			if (!isCompressed) {  return " P2SH: Only compressed key ";  }
			address[0] = 0x05;
			break;			
		case BECH32:
			if (!isCompressed) {	return " BECH32: Only compressed key ";	}
			hiiu_secp->GetHash160(type, isCompressed, pubKey, address + 1);

			char addr_bech32[128];
			segwit_addr_encode(addr_bech32, "bc", 0, address + 1, 20);
			return std::string(addr_bech32);
	}

	hiiu_secp->GetHash160(type, isCompressed, pubKey, address + 1);	

	sha256_checksum(address, 21, address + 21);

	std::string addr = EncodeBase58(address, address + 25);

	delete hiiu_secp;
	return addr;

}

//--------------------------------------------------------------------
std::string Bitcoin::hash160ToAddr(uint32_t* _hash160, bool isCompressed){
	// _hash160[] = 1784337440; 1882531190; -1293883124; -1810242794; -1514644173;
	// 	h[0]      = 1784337440; 1882531190;  3001084172;  2484724502;  2780323123 ;
	Secp256K1* hiiu_secp = new Secp256K1();   
	hiiu_secp->Init();	

	unsigned char address[25];
	address[0] = 0x00;
	memcpy(address + 1, _hash160, 20);
	sha256_checksum(address, 21, address + 21);
	std::string addr = EncodeBase58(address, address + 25);

	delete hiiu_secp;
	return addr;
	
}
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//--------------------------------------------------------------------


int main(){
	std::cout<<std::endl<<" -------------------------------------------- ";

	bool isCompressed = 1;
	Bitcoin bitcoin; 

	//priv --> Addr
	std::string addr = bitcoin.privToAddr(P2PKH, "2832ed74f2b5e35ee", isCompressed);
	std::cout << std::endl << "privToAddr addr : " << addr;
	addr = bitcoin.privToAddr(P2PKH, "2832ed74f2b5e35ee", !isCompressed);
	std::cout << std::endl << "privToAddr addr : " << addr;
	addr = bitcoin.privToAddr(P2SH, "2832ed74f2b5e35ee", isCompressed);
	std::cout << std::endl << "privToAddr addr : " << addr;
	addr = bitcoin.privToAddr(BECH32, "2832ed74f2b5e35ee", isCompressed);
	std::cout << std::endl << "privToAddr addr : " << addr;
	std::cout << "\n\n -------------------------------------------- \n";



	//priv --> hash160 
	uint32_t h[5]; 
	bitcoin.privToHash160(P2PKH, "2832ed74f2b5e35ee", h, isCompressed);
	bitcoin.privToHash160(P2PKH, "2832ed74f2b5e35ee", h, !isCompressed);
	bitcoin.privToHash160(P2SH, "2832ed74f2b5e35ee", h, isCompressed);
	bitcoin.privToHash160(BECH32, "2832ed74f2b5e35ee", h, isCompressed);
	std::cout<<"\n\n -------------------------------------------- \n";



	//priv --> Pubkey
	Point publicKey = bitcoin.privToPubkey("2832ed74f2b5e35ee");
			// std::cout<<std::endl<<"\n.\n.\n.\n---- public key : " << publicKey.toString();
	//pubkey --> Hash160
	uint32_t hash160[5]; 
	bitcoin.pubkeyToHash160(P2PKH, publicKey, hash160, isCompressed);
	bitcoin.pubkeyToHash160(P2PKH, publicKey, hash160, !isCompressed);
	bitcoin.pubkeyToHash160(P2SH, publicKey, hash160, isCompressed);
	bitcoin.pubkeyToHash160(BECH32, publicKey, hash160, isCompressed);			
			// for (int i = 0; i < 5; i++){	printf("\n --- _hash160[] : %d ", get_h[i]);
	std::cout<<"\n\n -------------------------------------------- \n";




	// pubkey --> addr 
	std::string addr1 = bitcoin.pubkeyToAddr(P2PKH, publicKey, isCompressed);
	std::cout << "\npubkeyToAddr addr1 : "<<addr1;
	addr1 = bitcoin.pubkeyToAddr(P2PKH, publicKey, !isCompressed);
	std::cout << "\npubkeyToAddr addr1 : "<<addr1;
	addr1 = bitcoin.pubkeyToAddr(P2SH, publicKey, isCompressed);
	std::cout << "\npubkeyToAddr addr1 : "<<addr1;
	addr1 = bitcoin.pubkeyToAddr(BECH32, publicKey, isCompressed);
	std::cout << "\npubkeyToAddr addr1 : "<<addr1;
	std::cout<<"\n\n -------------------------------------------- \n";


	printf("\n\n\n");
	return 0;
};
