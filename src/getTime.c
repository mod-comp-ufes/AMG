# include <sys/time.h>
# include "getTime.h"

double get_time () {
	struct timeval tv;
	gettimeofday(&tv, 0);
	return (double)(tv.tv_sec*100.0 + tv.tv_usec/10000.0);
};
