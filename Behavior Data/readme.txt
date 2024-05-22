Raw behavioral data for "The role of motor cortex in motor sequence execution depends on demands for flexibility"



Data is organized at matlab structures, where each entry is a session from a different training session, recorded over the entire course of the experiment

Description of relevant fields:
 -startTime: matlab datetime of when the session started
 -name: name of rat
 -pokeTimes: teensy time (in ms from start of session) when a lever was pressed. each element in cell array is a different consecutive trial
 -pokeNames: string of lever identity ('L','C','R') for detected lever presses. each element in cell array is a different consecutive trial
 -cuedNames: string of lever identity ('L','C','R') for when visual cue was presented.
 -Hit: vector of whether reward was delivered on a given trial
 -blocknumRepair: 0-5 vector indicating place in trial block.
 -protocol: integer ID of which session type the rat is in. 0-6 are training sessions. 7 is full CUE and WM trials. 8 is overtrained trials. N.b. some AUTO-only rats had a distinct set of protocols, and protocol 6 is overtrained trials
 -extraPokes: times of lever presses (in ms from start of session) performed during inter trial interval
 -extraPokesNames: identity of lever presses performed during inter trial interval

Session IDs of lesions

E1-Rat20 
 Unilateral: 702-743
 Bilateral: 837-873
J2-Rat17
 Unilateral: 1069-1119
 Bilateral: 1284-1319
J8-Rat7
 Unilateral: 810-844
 Bilateral: 960-995
F4-Rat43
 Unilateral: 657-700
 Bilateral: 805-830
E4-Rat39
 Unilateral: 1336-1368
 Bilateral: 1575-1607
L3-Rat34
 Unilateral: 1411-1452
 Bilateral: 1810-1841
J4-Rat126
 Unilateral: 1454-1478
 Bilateral: 1622-1646


L1-Rat52
 Unilateral: 300-324
 Bilateral: 391-420
E8-Rat82
 Unilateral: 349-376
 Bilateral: 457-496
J3-Rat88
 Unilateral: 281-310
 Bilateral: 405-428
L5-Rat95
 Unilateral: 270-293
 Bilateral: 394-418
J8-Rat79
 Unilateral: 559-590
 Bilateral: 957-991
L1-Rat77
 Unilateral: 709-741
 Bilateral: 862-898
