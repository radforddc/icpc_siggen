# Makefile for signal generation from PPC detectors
#   - uses .c library codes by Karin Lagergren, heavily modified by David Radford
#
# [-lreadline option required for readline, addhistory...]

CC = gcc
CPP = g++
CFLAGS = -O3 -Wall
RM = rm -f

# common files and headers
mk_signal_files = calc_signal.c cyl_point.c detector_geometry.c fields.c point.c read_config.c
mk_signal_headers = calc_signal.h cyl_point.h detector_geometry.h fields.h mjd_siggen.h point.h

All: stester mjd_fieldgen mass icpc_fieldgen_opt icpc_opt

extras: speed dtc csv2dat aoe_dist

# interactive interface for signal calculation code
stester: $(mk_signal_files) $(mk_signal_headers) signal_tester.c
	$(CC) $(CFLAGS) -o $@ $(mk_signal_files) signal_tester.c -lm -lreadline

mjd_fieldgen: mjd_fieldgen.c read_config.c detector_geometry.c mjd_siggen.h detector_geometry.h
	$(CC) $(CFLAGS) -o $@ mjd_fieldgen.c read_config.c detector_geometry.c -lm

icpc_fieldgen_opt: icpc_fieldgen_opt.c read_config.c detector_geometry.c mjd_siggen.h detector_geometry.h
	$(CC) $(CFLAGS) -o $@ icpc_fieldgen_opt.c read_config.c detector_geometry.c -lm

icpc_opt: icpc_opt.c
	$(CC) $(CFLAGS) -o $@ icpc_opt.c -lm

mass: mass.c read_config.c mjd_siggen.h
	$(CC) $(CFLAGS) -o $@ mass.c read_config.c -lm


speed: $(mk_signal_files) $(mk_signal_headers) speed.c
	$(CC) $(CFLAGS) -o $@ $(mk_signal_files) speed.c -lm

dtc: $(mk_signal_files) $(mk_signal_headers) dtc.c
	$(CC) $(CFLAGS) -o $@ $(mk_signal_files) dtc.c -lm

aoe_dist: $(mk_signal_files) $(mk_signal_headers) aoe_dist.c
	$(CC) $(CFLAGS) -o $@ $(mk_signal_files) aoe_dist.c -lm -lreadline

csv2dat: csv2dat.c
	$(CC) $(CFLAGS) -o $@ csv2dat.c


FORCE:

clean: 
	$(RM) *.o core* *[~%] *.trace
	$(RM) stester mjd_fieldgen mass dtc
