
/* timers.c */

long long readTSC(void);	/* not sure if we keep this public */
void init_timers(int n);
void stamp_timers(int i);
long long diff_timers(int i, int j);
