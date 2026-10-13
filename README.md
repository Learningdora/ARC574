 # ARC 574 Midterm: Finding the maximum lean angle of a corrugated wall

The project's objective is to find the maximum lean angle of a wall that can be created by offsetting the bricks of each course to some overhang values such that the stability of the whole wall is maintained. This considers how the robot is building the brick wall, stacking one brick on top of another. As such, the stability of the courses assemblage k+1 to N (with N being the last course) must be accessed in relation to course k.
The wall follows a running-bond configuration.
The objective of the project and the logic of the code is prompted by the designing team, and the code was implemented using Deepseek-flash-v4 (access provided by NRP). Claude Opus 5-5 was used to review the code and document version control.

## Requirements
- Rhino 8 + Grasshopper

## Scripts:
- curve_generate.py : Define base curve
- brick_corbel_column.py : Create individual columns with max lean angle
- brick_wall_populate.py : Append columns into wall assembly

## Units: 
- mm in GH/Rhino model

  






