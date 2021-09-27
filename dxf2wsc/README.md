## What's about ![Here should be the logo](doxygen/wsclogo.png)

### Background

The boiler water circulation program "wsc.exe" calculates the flow in the network of tubes that represent a boiler. It calculates the densities and resistances. The resulting flow will be compared to some criteria for safe flow. 

The calculation is done on a network of tubes limited only by the storage capacity of the computer. Test cases with 3000+ tubes have been calculated. 

Each tube has some properties like 3D position of start and end point, diameter, wall thickness, heat absorption, bend radius, orifice etc.

Setting up a data file for several 1000 tubes is quite time consuming. 

Therefore a kind of pre-processor is used. The name of the program is "dxf2wsc.exe". The whole input (as well as output) is based on geometry data in .dxf format. Any program that can generate 3D-line .dxf data and have the availability of a big number of layers with long layer names can be used. Unfortunately, most of open source design or drawing programs fail in one of the preconditions (AFAIK). 

The geometry has to be set up in the CAD program and each tube has to be on one layer with the additional data. Similar tubes can use the same layer. 

"dxf2wsc" reads the dxf-file and writes the data file for the circulation program. If there are problems with the dxf-file those problems will be given as error messages. If there is a problem with the geometry, like tube ends that are not connected to other tubes, a .dxf-file will be produced. It can be overlayed with the original drawing (model) and indicates the position of the problem. 

Finally, the general data like drum pressure, calculation method, tube roughness, convergence limit, maximum iterations etc. have to be given. 

In the second step the calculation program "wsc.exe" can be started. 

###The file structure

The program files "dxf2wsc" as well as "wsc" should be in same directory. The name or location of this main directory doesn't matter. Let's call it "Boiler Circulation". What matters is the directory structure below. 
We should have one subdirectory with the name "WSCDATA". Within this directory each project is situated in a different subdirectory. Lets call it "Boiler1", "Boiler2" (certainly you would prefer to use more descriptive name like the real project name. The directory "Boiler1" etc. contains all the data, like .dxf file of input data, the circulation input data file, the result files, and the criteria .dxf files (and some other intemediate files that helps tracking, if some problems occur).   

**file structure** 

|Directory  | Containing |
| ----------- | ----------- |
|Boiler Circulation|   WSCDATA, dxf2wsc, wsc, (maybe some dll, to fit for the compiler, f.i. MinGW) |
|WSCDATA | Boiler1, Boiler2, ... |
|Boiler1 | Boiler1.dxf, Later on: Boiler1.dat, Boiler1.res (the result txt file) and some .dxf files|
|Boiler2 | Boiler2.dxf, Later on: Boiler2.dat, Boiler2.res (the result txt file) and some .dxf files|
| ... | ... |

###The program

"dxf2wsc.exe" as well as the main calculation program "wsc.exe" are console applications. (Up to now) there is no window wrapper. 

To start the programs you can simply click on them in file explorer and start. A console will pop up and the program will guide you through the process. After finish the console will close. This means that any final messages are lost. My recommendation is opening a console (cmd.exe or PowerShell) and start the program. In Windows file explorer go to directory "Boiler Circulation". In the ribbon on top of this window click to "File" and "Open Windows PowerShell" -> "Open Windows PowerShell". Normally the PowerShell has a blue background.  At the command prompt ">" type "./dxf2wsc". This will start the program. By the unix-like command structure of PowerShell you are required to add "./". If you use the command shell (cmd.exe, black background) this is not needed. 

At first you will be asked for the project name (in our example "Boiler2"). This Project name can also be given as a parameter, i.e. you can start the program "./dxf2wsc Boiler2".

The program opens "\WSCDATA\Boiler2\Boiler2.dxf" and reads the data. 

If there is a problem with the data file, you will get an error message indicating approximately the line number in .dxf file. I never got it straight on the first trial. Additionally, there will be a file "Boiler2_Error.dxf". This file can be overlayed with the original drawing (Boiler2.dxf) and shows lines from steam drum center to the point where the error was found. The error should be corrected. This is easier than indentifying a line in original dxf file.  

The program looks for following errors: 
- zero length line or line shorter than 10 mm, 
- double lines (same starting and end points), 
- points, that only exist in one line,(each tube has to be connected to at least one other tube)
- intersection of tubes not at start or end point,
- wrong tube outside diameter or tube thickness (resulting in too small inside diameter) 

If there are errors the program stops and an indication of the position can be found in "Boiler2_Error.dxf".

For bidrum boilers see below.

If there is no error the basic data for the calculation are put in. 

1. maximum number of iterations
2. tolerance 
3. control parameter for protocol file
4. calculation method for 2-phase density and pressure drop
5. drum pressure in MPag 
6. tube roughness in mm
7. LevelFactor between highest and lowest heated point 
8. Water level in drum, difference to drum center in mm
9. Resistance of drum internals in kPa
10. Circulation ratio for start values

In detail:

1.maximum number of iterations\n
the program "wsc.exe" calculates the flow in the different tubes in an iteration. This gives the maximum number of iterations. Usually the result or close to the result is got within 100 iteration steps.   

2.tolerance\n
If the difference of sum of all flows between this and the last iteration step is below this tolerance, the calculation will be stopped. Another stopping criterion is the maximum percentual flow difference in the branches. A number between 0.1 and 1. is a good starting point.
 
3.control parameter for protocol file\n 
If there is a problem during calculation the program stops (leaving sometimes cryptic messages like "Error in H2O.enthalpy t < 0 degC") This example means that the function for water and steam properties should calculate the enthalpy at a temperature below freezing point. This function doesn't know which other function called it and why and where it got a temperature below 0 degC. For this reason part of the calculation process can be written to a protocol file, named f.i. "Boiler2.pro". As writing everything can produce **HUGE** files that are even unreadable, certain subset can be written.  
| show | index |
| ----------- | ----------- |
|nothing| 0 |
|info about mesh/network | 1 |
|info about calculation of enthalpy at branches inlet | 2 |
|info about pressure drop in tubes | 3 |
|info about pressure drop in branches | 5 |
|info about nodal pressure calculation | 6 |
|info and .dxf about flow for next iteration step | 7 |
|info about reversal of flow direction in branches and tubes | 8 |
|all info (careful, protocol file can get **HUGE**) | 10 |
|more details about mesh/network and path finding | 11 |
|more details about pressure drop in tubes | 13 |
|more info reversal of flow direction in branches and tubes | 18 |
|everything (careful, protocol file can get **REALLY HUGE**) | 20 |

4.calculation method for 2-phase density and pressure drop\n
Unfortunately there is no universal formula to calculate the density and dynamic pressure drop of 2-phase water-steam flow. Different authors came up with different formulae and consequently the results also differ. Several models were built in. Some authors only handle void fraction, density. In this case the dynamic pressure drop is calculated according VDI Heat Atlas (2nd edition).

| calculation method | index |
| ----------- | ----------- |
| Jirous/Jirous | J |
| Rouhani/Becker | R |
| Chexal/VDIHeatAtlas | E |
| Welzer/Welzer | W |
| Woldesemayat/VDIHeatAtlas | G |

5.drum pressure in MPag

6.tube roughness in mm\n
It's the tube roughness, new tubes with magnetite layer it's ~0.04 mm. For corroded tubes this can go up.\n 

7.LevelFactor between highest and lowest heated point, 0 ... 1\n 
In this program version downcomers are not heated (bidrum see below). As there are other connecting tubes that are also not heated there has to be another criterion for downcomer branches. They have to pass through a horizontal plane which elevation is given by this factor. Plane passing through lower distribution headers = 0. Plane passing through upper headers = 1. The factor has to be chosen that no other un-heated connection tubes pass through this plane.

8.Water level in drum, difference to drum center in mm\n
The water level in steam drum is given as difference to drum center. If water level is above to number is positive, if below negative.

9.Resistance of drum internals in kPa\n
Drum internals like baffles or cyclones cause a pressure drop that influences the water circulation.

10.Circulation ratio for start values \n
The starting values the flow in different branches are based on the circulation ratio. In this case an overall circulation ratio is used. A number between 10 and 30 would be a good initial guess. Certainly, the circulation ratio for each tube will be calculated for each tube during the iteration. 

After getting all input data the program writes the data to a data file in WSC data format, f.i. "Boiler2.dat". Each tube and each start/end point get a number that is further used. To identify the position of say "tube 537" or "point 486" another .dxf file "Boiler2_TubeNo.dxf" will be written. You can open this file with your CAD Program. The tube number is given in the layer of this tube (each tube has a separate layer). The points are represented by small crosses. Also the number is given in the layer.  

###Bidrum boiler (!experimental!)

"wsc.exe" can also calculate bidrum boilers. This is regarded as experimental because if the boiler has only heated downcomers the number of downcomers have to be given by the inlet enthalpy in the downcomer tubes (see below). This is not balanced with feed water flow and enthalpy. The number of downcomers is not calculated automatically.  

Each tube has to be connected to at least 2 other tubes (inlet (start point) and outlet (end point)). Even the drum(s) count as tubes or a number of tubes or better pseudo tubes. There should be a connection tube between the point where the tube enters the drum and the drum center. In bidrum boilers with 1000+ tubes between the drums this means 2000+ drum tubes. 

To ease the set up of the model of bidrum boiler in CAD program, those drum tubes can be omitted. It has only to made sure that the start/end points lay within the drum shells. 

If "dxf2wsc" detects a big number of loose ends, it asks if we have a bidrum boiler. If yes, it ask for 

1. steam drum outside diameter in mm
2. steam drum shell thickness in mm
3. absolute x-coordinate of mud drum center point in mm
4. absolute y-coordinate of mud drum center point in mm
5. absolute z-coordinate of mud drum center point in mm
6. mud drum outside diameter in mm
7. mud drum shell thickness in mm

The drum pseudo tubes will be generated.

##The .dxf file

The .dxf format was introduced to exchange drawing data between CAD programs and between CAD and other programs. It is a readable text format with a big number of lines.  

As the CAD program can only handle geometric objects like lines and arc etc. the geometry, i.e. start and end points of lines(tubes) /arcs(bends, elbows) additional data relating to a tube has to be given elsewhere or in a different way. As each line/arc is related to a layer, providing this information by the layer is an easy way. Certainly, lines can have attributes but this is more complicated. 

The boiler is modelled as a 3D wire model. Basically, it is the centerlines of all tubes. It has to be set up in 3D. If there is only a 2D drawing of the boiler body, this need to be converted to 3D model. If the boiler is set up as a 3D solid model, the tube centerlines need to be extracted. What is needed is a network of lines/arcs in 3D space. All dimensions in mm. 

The colors and line types are ignored. 

The model should be saved in 2013 dxf format. Other formats are not tested by now. 

The z-coordinate should represent the elevation. 

The steam drum centerpoint should be at x = 100 000 mm, y = 100 000 mm, z = 100 000 mm. Just to be recognized as steam drum. The steam drum has some diamter and lenght. It should be represented by it's centerpoint. The tubes/pipes are connected to drum at points away from centerpoint. There should be connection pseudo tubes between tube/pipe start/end point to drum center. Those pseudo tubes get diameter and thickness from drum. For bidrum boilers those pseudo tubes can be generated by the program (see above). 

As mentioned the layer name contains the further information for the tubes/pipes. It organised as follows: 

"tube name_outside diameter_wall thickness_heat absorption per m_number of identical tubes_additional resistance factor_diameter of orifice at inlet_diameter of orifice at outlet_factor for position of heating_tube inlet enthalpy"

The different numbers are separated by underscores. "front wall_60.3_4.5_36_1_0_0_0_1_0" or in short "front wall_60.3_4.5_36" is valid. "_1_0_0_0_1_0" are standard value and can be omitted. If there is a number that doesn't contain a standard value the preceeding numbers have to be given. "downcomer_50.8_3.2_5_1_0_0_0_1_-2" can not be shortend. 

| number | description | standard value| mandatory|
| ----------- | ----------- |----------- |----------- |
|tube name|any name helping to identify the tube | none | yes |
|outside diameter| tube outside diameter in mm | none | yes |
|wall thickness| tube wall thickness in mm | none | yes |
|heat absorption per m | in kW see note | 0 | no |
|number of identical tubes | see note | 1 | no |
|additional resistance factor | additional resistance factor | 0 | no |
|diameter of orifice at inlet | diameter of orifice at inlet in mm | 0 | no |
|diameter of orifice at outlet | diameter of orifice at outlet in mm | 0 | no |
|factor for position of heating | see note | 1 | no |
|tube inlet enthalpy | in kJ/kg see note | 0 | no|

*heat absorption per m*: in some areas the heat absorption per m2, f.i. tube walls, is constant but there are tubes of different length. To make it easier the heat absoption per m is given and later multiplied with the length by the program. This helps to keep the number of different layers low.

*number of identical tubes* Basically for the model several real tubes can be replaced by one "model" tube. In this case this number is the number of tubes. "heat absorption per m" should also relate to the real numbers. Example: 2 real tubes are replaced by 1 tube. this number is 2. if each real tube has a heat load of 4 KW/m the heat absortion per m should be 8 kW/m. 

*additional resistance factor* dimensionless number that multiplied with velocity head give a pressure drop.

*diameter of orifice at in/outlet* hydraulic diameter of orifice or smaller hole in connecting header.

*factor for position of heating* the total density of steam-water mixture in a non-uniformly heated tube depends on the position of the heating. This is not taken into account but kept for future program versions.

*tube inlet enthalpy* this is for bidrum boilers. To avoid steaming in heated downcomers we have to make sure that at no place the water saturation enthalpy is exceeded. As by the tube length downwards the pressure goes up this is a tricky balance. To help in this matter the feed water line is close to the inlet of the downcomer tubes resulting in an inlet enthalpy below drum water saturation enthalpy.   
This number is the difference to drum water enthalpy. Bidrum heated downcomer should have a negative number, say -5 or -10. This number need to be verified after calculation of circulation.  

###The data file

"wsc.exe" needs its own input format. "dxf2wsc.exe" reads the .dxf file and finally writes the file, f.i."Boiler2.dat".

This file is a text file that can be opened with each text editor. If only "number of iterations" should be changed it's a hazzle to go through all the input data in "dzf2wsc" again just to change this particular number. Certainly, if geometry data should be changed it is better to start all over in the CAD program. 

**the layout of data file**

First line contains the project name incl. paths.

Second line 

1. maximum number of iterations
2. tolerance 
3. control parameter for protocol file
4. calculation method for 2-phase density and pressure drop

Third line 

1. drum pressure in MPag 
2. tube roughness in mm
3. LevelFactor between highest and lowest heated point 
4. Water level in drum, difference to drum center in mm
5. Resistance of drum internals in kPa
6. Circulation ratio for start values

following lines 

| number of point | x-coordinate | y-coordinate | z-coordinate | 

"TubeData" is the separator between point and tube data

each tube data consist of 2 lines

first line: tube name

second line 

| number | description | 
| ----------- | ----------- |
| PointIn | number of inlet point |
| PointOut | number of outlet point |
| dOut | tube outside diameter |
| Thickness | tube wall thickness
| RadiusBend | bending radius |
| NoParallel | number of identical (parallel) tubes |
| ksiAdd | additional resistance factor |
| q | heat absorption of tube in kW | 
| FactHeat | factor for heating position |
| DiaOrificeIn | diameter of orifice at inlet of tube |
| DiaOrificeOut | diameter of orifice at outlet of tube |
| EnthGiven | inlet enthalpy difference to downcomers in bidrum boilers |
| OCSStartAngle | start angle of arc |
| OCSEndAngle | end angle of arc |
| OCSCenterX | x-coordinate center of arc |
| OCSCenterY | y-coordinate center of arc |
| OCSCenterZ | z-coordinate center of arc |
| Nx | x-component of normal vector |
| Ny | y-component of normal vector |
| Nz | z-component of normal vector |

OSC... to Nz are data needed to write dxf files for the final results. They are related to arc which have an object coordinate system. Those numbers are available from the original dxf file and are transferred not to be calculated again. 

