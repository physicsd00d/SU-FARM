This folder contains a few sample compact envelopes.  All envelopes are generated using a 5-minute reaction time and instantaneous space-vehicle health updating.

LynxMII_OK: This is a Lynx-like vehicle (HTHL) operating out of the Oklahoma Spaceport.  Shape boundary corresponds to a cumulative probability of debris strike constrained to be less than 1e-7.

Falcon9_Cape_Thresh: This is Falcon9-like launch to the ISS from Cape Canaveral.  Cumulative probability of debris strike < 1e-7.

Falcon9_Cape_CondSafe: This is Falcon9-like launch to the ISS from Cape Canaveral.  DOES NOT USE a probabilistic safety threshold.  Instead, this restricts the areas of the airspace where the debris will be found within 5 minutes (the reaction time) for ANY time of failure.  In other words, if aircraft avoid this shape, their risk is effectively zero.

* Each subsection has a 6-line preamble
NAMEX_Y                               # This is the NAME of the SUA, X is the timestep (starts at zero), and Y is the index for shapes within timestep X
                                      # Empty space
0                                     # I'm not sure.  Make this zero.
0 99900 54600 56580                   # LowerAlt UpperAlt StartTime EndTime, units are ft. and seconds-since-midight
0 0 0 1                               # I'm not sure.
11                                    # Number of lat/lon points which will follow immediately

The preamble section looks better when you download the file.  I dunno what's wrong with github not displaying the line breaks or spaces.  Good luck!
