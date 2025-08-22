KPL/FK
 
   FILE: lander_guess.tf
 
   This file was created by PINPOINT.
 
   PINPOINT Version 3.3.0 --- December 13, 2021
   PINPOINT RUN DATE/TIME:    2024-12-26T12:00:01
   PINPOINT DEFINITIONS FILE: lander.def
   PINPOINT PCK FILE:         pck00010.tpc
   PINPOINT SPK FILE:         lander_guess.bsp
 
   The input definitions file is appended to this
   file as a comment block.
 
 
   Body-name mapping follows:
 
\begindata
 
   NAIF_BODY_NAME                      += 'LANDER'
   NAIF_BODY_CODE                      += 301001
 
\begintext
 
 
   Reference frame specifications follow:
 
 
   Topocentric frame LANDER_TOPO
 
      The Z axis of this frame points toward the zenith.
      The X axis of this frame points North.
 
      Topocentric frame LANDER_TOPO is centered at the
      site LANDER, which has Cartesian coordinates
 
         X (km):                  0.3489173019355E+03
         Y (km):                  0.9349210927886E+02
         Z (km):                  0.1699433641515E+04
 
      and planetodetic coordinates
 
         Longitude (deg):        15.0000000000000
         Latitude  (deg):        78.0000000000000
         Altitude   (km):         0.0000000000000E+00
 
      These planetodetic coordinates are expressed relative to
      a reference spheroid having the dimensions
 
         Equatorial radius (km):  1.7374000000000E+03
         Polar radius      (km):  1.7374000000000E+03
 
      All of the above coordinates are relative to the frame IAU_MOON.
 
 
\begindata
 
   FRAME_LANDER_TOPO                   =  1301001
   FRAME_1301001_NAME                  =  'LANDER_TOPO'
   FRAME_1301001_CLASS                 =  4
   FRAME_1301001_CLASS_ID              =  1301001
   FRAME_1301001_CENTER                =  301001
 
   OBJECT_301001_FRAME                 =  'LANDER_TOPO'
 
   TKFRAME_1301001_RELATIVE            =  'IAU_MOON'
   TKFRAME_1301001_SPEC                =  'ANGLES'
   TKFRAME_1301001_UNITS               =  'DEGREES'
   TKFRAME_1301001_AXES                =  ( 3, 2, 3 )
   TKFRAME_1301001_ANGLES              =  (  -15.0000000000000,
                                             -12.0000000000000,
                                             180.0000000000000 )
 
\begintext
 
 
Definitions file lander.def
--------------------------------------------------------------------------------
 
begindata
 
    SITES =  'LANDER'
 
    LANDER_CENTER  = 301
    LANDER_FRAME   = 'IAU_MOON'
 
    LANDER_IDCODE  = 301001
    LANDER_LATLON  = (78, 15, 0)
    LANDER_UP      = 'Z'
    LANDER_NORTH   = 'X'
 
 
 
 
begintext
 
[End of definitions file]
 
