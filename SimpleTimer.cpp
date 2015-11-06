#include "SimpleTimer.h"

FLOAT timeDifference(const struct timeval & t1, const struct timeval & t2){
    long seconds, useconds;
    seconds = t2.tv_sec - t1.tv_sec;
    useconds = t2.tv_usec - t1.tv_usec;
    return (FLOAT)seconds + useconds / 1.0e6;
}

SimpleTimer::SimpleTimer(){}

void SimpleTimer::start(){
    gettimeofday(& _startingTime, NULL);
}

FLOAT SimpleTimer::getTimeAndRestart(){
    gettimeofday(& _currentTime, NULL);
    FLOAT returnValue = timeDifference(_startingTime, _currentTime);
    _startingTime = _currentTime;
    return returnValue;
}

FLOAT SimpleTimer::getTimeAndContinue(){
    gettimeofday(& _currentTime, NULL);
    return timeDifference(_startingTime, _currentTime);
}
