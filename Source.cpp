#include "C:/Users/Casper/Desktop/C++/SunMoon.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <limits>
#include <string>
#include <iostream>
#include <time.h>

using namespace Astronomy;
static double m_MoonAz;
static double m_MoonAlt;

int main(void) {
	struct tm m_datumzeit;
	int m_Zone;
	double latitude = 55.782604799999994;
	double longitude = 12.528025600000001;
	double deltaT = 65.0;
	while (1) {
		time_t now = time(NULL);
#ifdef _MSC_VER
		struct tm gmt;
		memset(&gmt, 0, sizeof(struct tm));
		memset(&m_datumzeit, 0, sizeof(struct tm));
		localtime_s(&m_datumzeit, &now);
		gmtime_s(&gmt, &now);
#else
		m_datumzeit = *localtime(&now);
		struct tm gmt = *gmtime(&now);
#endif
		time_t utc = mktime(&gmt);
		if (gmt.tm_isdst == 1) utc -= 3600;
		m_Zone = (int)(difftime(now, utc) / 3600.0);
		Astronomy::SunMoon::SunMoon(latitude,longitude, &m_datumzeit, deltaT, m_Zone);

		&Astronomy::SunMoon::GetMoonAlt;
		&Astronomy::SunMoon::GetMoonAz;
		//std::cout << "MoonAlt:  "<< m_MoonAlt << "		Moonz:  " << m_MoonAz << "\n";
	}

	return 100;
}