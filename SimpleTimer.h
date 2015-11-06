#ifndef _SIMPLE_TIMER_H_
#define _SIMPLE_TIMER_H_

#include <sys/time.h>
#include "Definitions.h"


/** Timer. Used for scaling measurements.
 */
class SimpleTimer {

    public:
        SimpleTimer();

        /** Starts the timer */
        void start();

        /** Gets the elapsed time since last restart, and restarts the measurement
         * @return Elapsed time since last restart
         */
        FLOAT getTimeAndRestart();

        /** Gets the elapsed time since last restart.
         * @return Elapsed time since last restart
         */
        FLOAT getTimeAndContinue();

    private:

        struct timeval _startingTime;
        struct timeval _currentTime;

};

#endif
