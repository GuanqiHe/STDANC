#pragma once

#include "msgpack.hpp"

typedef std::vector<double> array_t;

struct datapack_t
{
	array_t t;
	array_t y;
	array_t d;
	array_t u;
	array_t phi;
	MSGPACK_DEFINE_MAP(t, y, d, u);
};

extern datapack_t data;