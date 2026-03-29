#ifndef _OPERATIONS128_H_
#define _OPERATIONS128_H_

#include <stdint.h>

typedef struct {
	uint64_t low;
	uint64_t high;
} uint128;

inline int _isequal_128(uint128 v1, uint128 v2) {
	if (v1.low==v2.low && v1.high==v2.high) {
		return 1;
	} else {
		return 0;
	}
}

inline uint128 _shift_left_128(uint128 val, char n) {
	uint128 ret;
	if (n==0) {
		ret.high = val.high;
		ret.low = val.low;
	} else if (n>0 && n<=64) {
		ret.high = (val.high << n) | (val.low >> (64-n));
		ret.low = val.low << n;
	} else if (n>64 && n<128) {
		ret.high = val.low << (n-64);
		ret.low = 0;
	} else {
		ret.high = 0;
		ret.low = 0;
	}
	return ret;
} 

inline uint128 _shift_right_128(uint128 val, char n) {
	uint128 ret;
	if (n==0) {
		ret.high = val.high;
		ret.low = val.low;
	} else if (n>0 && n<=64) {
		ret.high = val.high >> n;
		ret.low = (val.low >> n) | (val.high << (64-n));
	} else if (n>64 && n<128) {
		ret.high = 0;
		ret.low = val.high >> (n-64);
	} else {
		ret.high = 0;
		ret.low = 0;
	}
	return ret;
} 

inline uint128 _and_128(uint128 v1, uint128 v2) {
	uint128 ret;
	ret.high = v1.high & v2.high;
	ret.low = v1.low & v2.low;
	return ret;
}

inline uint128 _or_128(uint128 v1, uint128 v2) {
	uint128 ret;
	ret.high = v1.high | v2.high;
	ret.low = v1.low | v2.low;
	return ret;
}

inline uint128 _not_128(uint128 v) {
	uint128 ret;
	ret.high = ~v.high;
	ret.low = ~v.low;
	return ret;
}

#endif