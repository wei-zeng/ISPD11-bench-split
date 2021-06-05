#ifndef _VIA_H_
#define _VIA_H_

#include "baseDB.h"

class Via {
public:
	Via() :
			_cap(0), _demand(0), _rippedUp(false) {
	}

	//  write
	void setCap(int area) {
		_cap = area;
	}
	void setLoc(int x, int y, int z) {
		_x = x;
		_y = y;
		_z = z;
	}
	void addDemand(int area) {
		_demand += area;
	}
	void ripUp() {
		_rippedUp = true;
	}
	void restore() {
		_rippedUp = false;
	}
	void markOriginal() {
		_original = true;
	}
	void demarkOriginal() {
		_original = false;
	}
	//  read
	int cap() const {
		return _cap;
	}
	int demand() const {
		return _demand;
	}
	int of() const {
		return max(0, _demand - _cap);
	}
	Gcell getLoc() const {
		return Gcell(_x, _y, _z);
	}
	bool isRippedUp() const {
		return _rippedUp;
	}
	bool isOriginal() const {
		return _original;
	}
	bool operator==(const Via &other) const {
		return _cap == other._cap && _demand == other._demand;
	}
private:
	int _cap, _demand;
	int _x, _y, _z;
	bool _rippedUp, _original;
};

#endif
