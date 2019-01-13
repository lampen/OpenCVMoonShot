/*
Entfernen Sie folgende Informationen auf keinen Fall: / Do not remove following text:
Source code based on the javascript by Arnold Barmettler, www.astronomie.info / www.CalSky.com
based on algorithms by Peter Duffett-Smith's great and easy book
'Practical Astronomy with your Calculator'.
*/
#pragma once
#include <time.h>
#include <iostream>
#include <map>

namespace Astronomy
{
	class TimeSpan
	{
		friend class SunMoon;
	private:
		int m_hh;
		int m_mm;
		int m_ss;
		double m_dv;

		TimeSpan(double tdiff);

	public:
		TimeSpan();
		~TimeSpan();

		int GetHour();
		int GetMinute();
		int GetSecond();

		float GetTotalHour();
		float GetTotalMinute();
		float GetTotalSecond();

		double GetFloat();
		operator double();

		friend std::ostream& operator<<(std::ostream& os, const TimeSpan& ts);
	};

	enum SIGN
	{
		SIGN_ARIES,			//!< Widder
		SIGN_TAURUS,		//!< Stier
		SIGN_GEMINI,		//!< Zwillinge
		SIGN_CANCER,		//!< Krebs
		SIGN_LEO,			//!< Löwe
		SIGN_VIRGO,			//!< Jungfrau
		SIGN_LIBRA,			//!< Waage
		SIGN_SCORPIO,		//!< Skorpion
		SIGN_SAGITTARIUS,	//!< Schütze
		SIGN_CAPRICORNUS,	//!< Steinbock
		SIGN_AQUARIUS,		//!< Wassermann
		SIGN_PISCES			//!< Fische
	};

	enum LUNARPHASE
	{
		LP_NEW_MOON,                //!< Neumond
		LP_WAXING_CRESCENT_MOON,    //!< Zunehmende Sichel
		LP_FIRST_QUARTER_MOON,      //!< Erstes Viertel
		LP_WAXING_GIBBOUS_MOON,     //!< Zunehmender Mond
		LP_FULL_MOON,               //!< Vollmond
		LP_WANING_GIBBOUS_MOON,     //!< Abnehmender Mond
		LP_LAST_QUARTER_MOON,       //!< Letztes Viertel
		LP_WANING_CRESCENT_MOON,    //!< Abnehmende Sichel
	};

	class SunMoon
	{
	private:
		double m_Lat;
		double m_Lon;
		struct tm m_datumzeit;
		int m_Zone;
		double m_DeltaT;
		double m_JD;
		TimeSpan m_GMST;
		TimeSpan m_LMST;
		double m_SunDistance;
		double m_SunDistanceObserver;
		double m_SunLon;
		TimeSpan m_SunRA;
		double m_SunDec;
		double m_SunAz;
		double m_SunAlt;
		double m_SunDiameter;
		TimeSpan m_SunAstronomicalTwilightMorning;
		TimeSpan m_SunNauticalTwilightMorning;
		TimeSpan m_SunCivilTwilightMorning;
		TimeSpan m_SunRise;
		TimeSpan m_SunTransit;
		TimeSpan m_SunSet;
		TimeSpan m_SunCivilTwilightEvening;
		TimeSpan m_SunNauticalTwilightEvening;
		TimeSpan m_SunAstronomicalTwilightEvening;
		SIGN m_SunSign;
		double m_MoonDistance;
		double m_MoonDistanceObserver;
		double m_MoonLon;
		double m_MoonLat;
		TimeSpan m_MoonRA;
		double m_MoonDec;
		double m_MoonAz;
		double m_MoonAlt;
		double m_MoonDiameter;
		TimeSpan m_MoonRise;
		TimeSpan m_MoonTransit;
		TimeSpan m_MoonSet;
		double m_MoonPhaseNumber;
		double m_MoonAge;
		LUNARPHASE m_MoonPhase;
		SIGN m_MoonSign;

	public:
		SunMoon();
		SunMoon(double latitude, double longitude, struct tm* datetime, double deltaT);
		SunMoon(double latitude, double longitude, struct tm* datetime, double deltaT, int zone);
		~SunMoon();

		void SetInput(double latitude, double longitude, struct tm* datetime, double deltaT);
		void SetInput(double latitude, double longitude, struct tm* datetime, double deltaT, int zone);

		/** Östl. geografische Länge [Grad] */
		double GetLat();

		/** Geografische Breite [Grad] */
		double GetLon();

		/** Datum & Uhrzeit */
		struct tm GetDateTime();

		/** Zeitdifferenz zu Weltzeit [h] */
		double GetZone();

		/** deltaT [sek] */
		double GetDeltaT();

		/** Julianisches Datum [Tage] */
		double GetJD();

		/** Greenwich Sternzeit GMST [h] */
		TimeSpan GetGMST();

		/** Lokale Sternzeit LMST [h] */
		TimeSpan GetLMST();

		/** Entfernung der Sonne (Erdmittelpunkt) [km] */
		double GetSunDistance();

		/** Entfernung der Sonne (vom Beobachter) [km] */
		double GetSunDistanceObserver();

		/** Eklipt. Länge der Sonne [Grad] */
		double GetSunLon();

		/** Rektaszension der Sonne [h] */
		TimeSpan GetSunRA();

		/** Deklination der Sonne [Grad] */
		double GetSunDec();

		/** Azimut der Sonne [Grad] */
		double GetSunAz();

		/** Höhe der Sonne über Horizont [Grad] */
		double GetSunAlt();

		/** Durchmesser der Sonne ['] */
		double GetSunDiameter();

		/** Astronomische Morgendämmerung [h] */
		TimeSpan GetSunAstronomicalTwilightMorning();

		/** Nautische Morgendämmerung [h] */
		TimeSpan GetSunNauticalTwilightMorning();

		/** Bürgerliche Morgendämmerung [h] */
		TimeSpan GetSunCivilTwilightMorning();

		/** Sonnenaufgang [h] */
		TimeSpan GetSunRise();

		/** Sonnenkulmination [h] */
		TimeSpan GetSunTransit();

		/** Sonnenuntergang [h] */
		TimeSpan GetSunSet();

		/** Bürgerliche Abenddämmerung [h] */
		TimeSpan GetSunCivilTwilightEvening();

		/** Nautische Abenddämmerung [h] */
		TimeSpan GetSunNauticalTwilightEvening();

		/** Astronomische Abenddämmerung [h] */
		TimeSpan GetSunAstronomicalTwilightEvening();

		/** Tierkreiszeichen */
		SIGN GetSunSign();

		/** Entfernung des Mondes (Erdmittelpunkt) [km] */
		double GetMoonDistance();

		/** Entfernung des Mondes (vom Beobachter) [km] */
		double GetMoonDistanceObserver();

		/** Eklipt. Länge des Mondes [Grad] */
		double GetMoonLon();

		/** Eklipt. Breite des Mondes [Grad] */
		double GetMoonLat();

		/** Rektaszension des Mondes [h] */
		TimeSpan GetMoonRA();

		/** Deklination des Mondes [Grad] */
		double GetMoonDec();

		/** Azimut des Mondes [Grad] */
		double GetMoonAz();

		/** Höhe des Mondes über Horizont [Grad] */
		double GetMoonAlt();

		/** Durchmesser des Mondes ['] */
		double GetMoonDiameter();

		/** Mondaufgang [h] */
		TimeSpan GetMoonRise();

		/** Mondkulmination [h] */
		TimeSpan GetMoonTransit();

		/** Monduntergang [h] */
		TimeSpan GetMoonSet();

		/** Mondphase */
		double GetMoonPhaseNumber();

		/** Mondalter [Grad] */
		double GetMoonAge();

		/** Mondphase */
		LUNARPHASE GetMoonPhase();

		/** Mondzeichen */
		SIGN GetMoonSign();

		friend std::ostream& operator<<(std::ostream& os, const SunMoon& sm);
	};

}
