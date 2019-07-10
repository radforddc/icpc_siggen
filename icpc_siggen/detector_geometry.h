/* detector_geometry.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module keeps track of the detector geometry
 */
#ifndef _DETECTOR_GEOMETRY_H
#define _DETECTOR_GEOMETRY_H

#include "mjd_siggen.h"
#include "point.h"
#include "cyl_point.h"

/* ouside_detector
   returns 1 if pt is outside the detector, 0 if inside detector
*/
int outside_detector(point pt, MJD_Siggen_Setup *setup);
int outside_detector_cyl(cyl_pt pt, MJD_Siggen_Setup *setup);

#endif /*#ifndef _DETECTOR_GEOMETRY_H*/
