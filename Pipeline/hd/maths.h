
#pragma once

template<typename T>
inline T div_round_up(T a, T b) { return (a-(T)1)/b+(T)1; }

template<typename T>
inline T min(T a, T b) { return a<b ? a : b; }
template<typename T>
inline T max(T a, T b) { return a<b ? b : a; }
