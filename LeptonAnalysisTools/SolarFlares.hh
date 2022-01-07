#ifndef SolarFlares_hh
#define SolarFlares_hh

bool IsTimeDuringSolarFlare(double timeStamp, double magnitudeTreshold = 1000.0, double excludeSecondsBefore = 2.0 * 3600.0, double excludeSecondsAfter = 22.0 * 3600.0);

#endif
