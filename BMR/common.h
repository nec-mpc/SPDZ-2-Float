// (C) 2018 University of Bristol, Bar-Ilan University. See License.txt

/*
 * common.h
 *
 */

#ifndef CIRCUIT_INC_COMMON_H_
#define CIRCUIT_INC_COMMON_H_

#include <string>

typedef unsigned long wire_id_t;
typedef unsigned long gate_id_t;

class Function {
	bool rep[4];
	int shift(int i) { return 4 * (3 - i); }
public:
	Function() { memset(rep, 0, sizeof(rep)); }
	Function(std::string& func)
	{
		for (int i = 0; i < 4; i++)
			if (func[i] != '0')
				rep[i] = 1;
			else
				rep[i] = 0;
	}
	Function(int int_rep)
	{
		for (int i = 0; i < 4; i++)
			rep[i] = (int_rep << shift(i)) & 1;
	}
	uint8_t operator[](int i) { return rep[i]; }
};

template <class T>
class CheckVector : public vector<T>
{
public:
	CheckVector() : vector<T>() {}
	CheckVector(size_t size) : vector<T>(size) {}
	CheckVector(size_t size, const T& def) : vector<T>(size, def) {}
#ifdef CHECK_SIZE
	T& operator[](size_t i) { return this->at(i); }
	const T& operator[](size_t i) const { return this->at(i); }
#else
	T& at(size_t i) { return (*this)[i]; }
	const T& at(size_t i) const { return (*this)[i]; }
#endif
};

#endif /* CIRCUIT_INC_COMMON_H_ */
