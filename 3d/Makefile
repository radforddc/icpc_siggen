# Makefile for signal generation from PPC detectors
#   - uses .c library codes by Karin Lagergren, heavily modified by David Radford
#
# [-lreadline option required for readline, addhistory...]

CC = gcc
CPP = g++
CFLAGS = -O3 -Wall
RM = rm -f

# common files and headers
mk_signal_files = field_init.c calc_signal.c cyl_point.c detector_geometry.c fields.c point.c \
                  read_config.c signal_calc_util.c
mk_signal_headers = calc_signal.h cyl_point.h detector_geometry.h fields.h mjd_siggen.h \
                    point.h signal_calc_util.h

All: stester mjd_fieldgen3d
# mass drift_time_map

# interactive interface for signal calculation code
stester: $(mk_signal_files) $(mk_signal_headers) signal_tester.c
	$(CC) $(CFLAGS) -o $@ $(mk_signal_files) signal_tester.c -lm -lreadline

mjd_fieldgen3d: fieldgen3d.c field_init.c point.c read_config.c signal_calc_util.c mjd_siggen.h point.h
	$(CC) $(CFLAGS) -o $@ fieldgen3d.c field_init.c point.c read_config.c signal_calc_util.c -lm

FORCE:

clean: 
	$(RM) *.o core* *[~%] *.trace
	$(RM) stester mjd_fieldgen3d
