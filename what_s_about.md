## What's about ![logo](WSCLogo.png)

### Background
This program is about steam boilers. I mean real boilers where you convert water to steam by using heat of any source.

Steam technology is about 300 years old but is still needed in some area of application. It is nearly impossible to convert say wood or garbage into electricity with a computer or other electronic gadgets alone. 

Beside more modern methods like solar panels, wind turbines, hydropower or tidal wave power stations there is still a demand for electricity from solid waste resources (like agricultural waste). Beside that, steam is even needed as a heating source in different industrial applications like refineries or food industry (and using electricity to produce the steam is far too expensive, just look at your electricity bill).

Nowadays, there are different boiler types best suited for the different applications. For low pressure (up to about 2MPa) fire-tube boilers (<a href=" https://en.wikipedia.org/wiki/Fire-tube_boiler">see Wikipedia</a>) are a budget solution for the steam demand. The heat transfer tubes have the hot gasses at the inside. On the outside is boiling water. They can be fired with liquid or gas fuel. In special cases designs to burn coal also exist. Above the low pressure range a different design is used where the water or steam is inside the tubes and the heat source (hot gas) is outside. Those are called water-tube boilers (<a href="https://en.wikipedia.org/wiki/Water-tube_boiler">see Wikipedia</a>) 

Again the application asks for different designs. Ranging from low pressure but high capacity boiler, f.i. in sugar factories to power stations with supercritical pressure and high temperatures to produce electricity in a very efficient way. 

Different types of water-tube boilers exist. A very popular type in certain areas of the world is the by-drum boiler. As the name suggests there are 2 drums that are connected by tubes that take up the heat from hot combustion gases. 

In other areas of the world a fully welded design with a single drum is more popular. The heated tubes are connected to bigger pipes and form a kind of box. Typically there is a set of pipes at the bottom of the mostly vertical tubes and a set at the top. The bigger pipes are called headers. 

As the heated tubes have to convert water to steam there has to be a supply of water to the tubes as well as a discharge of steam. In drum boilers the steam leaving the heated tubes contains still a lot of water. Actually there is 5% (some typical number, mass related) of steam in the water. This mixture goes to the drum where the steam will be separated and leaves the boiler. The water will be mixed with feed water that replaces the generated steam. Via mostly big sized pipes the water goes down to the lower headers and feed to the heated tubes. The pipes going down from drum to the lower part of boiler are called downcomers. They can be heated or unheated. Unheated downcomers should receive water at saturation temperature of the drum. Whereas heated downcomers should get water at lower temperature by mixing with lower temperature feed water. The lower temperature has to make sure that no significant amount of steam is generated in the downcomer tubes. A significant amount of steam will cause the circulation to break down at least in this particular downcomer tube/pipe.  

### The problem
Water at atmospheric pressure boils at 100 degC. To get steam at higher temperature (and containing more energy) the pressure has to be higher than atmospheric pressure. In the industrial range we have to deal with 2.0 to 12.0 MPa. (1 MPa= mega Pascal is about 10 times the atmospheric pressure) So, the tubes in water-tube boilers are pressurized from the inside. The tube material has to withstand this pressure. Each boiler undergoes a certification according the applicable codes and has to pass pressure test to assure safe operation. There are also safety devices like safety valves that relieve the pressure if the pressure goes over the operation limit. The safety valve(s) keep the boiler safe.... in theory. 

Tubes and pipes are made of steel. If you heat the steel it becomes weaker. At 900+ deg C the steel get so soft that you can forge it. The flames and combustion gases have temperatures well above. At outside of the tubes you can have high temperature. The keep the strength of the steel the overall temperature in the tube wall should be low. So, at the inside of the tube the temperature should be low enough to get a low overall temperature. You also need a high heat transfer at the inside of the tube. Boiling water or steam at a high velocity can assure a high heat transfer. In this case the outside metal temperature is only ~20 K (Kelvin) above the inside temperature and we are safe. 

The problem can be phrased quite simple: **Make sure that each tube contains enough water that no part of the tube is "dry"**. 

As the tubes are heated, water without replacement will evaporate and leaves the tube "dry". The above phrase should be qualified: Make sure that there is enough **water flow** that no part of the tube is "dry".

In vertical tubes the flow has to be high enough that a critical steam mass content (or steam quality) will not occur. In horizontal tubes it depends less on steam mass content but more on the question whether there is a separation of water and steam in the tube. The steam will stay at the upper part of the tube preventing sufficient cooling. However, in both cases a problem occurs if the flow is too low. 

A water-tube boilers consists of a number of connected tubes with individual diameters, lengths and resistances.  

One way to get enough flow is using a circulation pump on the downcomer part. This pump insures that there is enough flow in all tubes ....in theory.

As the tubes have different resistances the flow will distribute that the water has the easiest way of less resistance. Tubes with less resistance will get more flow than tubes with high resistance. 

A pump alone is no guarantee of trouble-free operation. You still have to make sure that there is enough flow in the high resistance tubes, f.i. by adding some resistance to the other tubes. 

As a pump is a part that is costly and needs electricity and maintenance, other ways to get sufficient circulation in the boiler are preferred. 

The downcomers are filled with water. Water has a certain density. On the other hand the heated tubes contain a mixture of water and steam. It's density (typically) is far below the water density. Connecting 2 vertical tubes or pipes result in a driving force for the circulation like a pump. The driving force is balanced by the resistance force. As no extra pump is needed this is called "natural circulation".

### The program
The program calculates the flow in the network of tubes that represent a boiler. It calculates the densities and resistances. The resulting flow will be compared to some criteria for safe flow. 

The calculation is done on a network of tubes limited only by the storage capacity of the computer. Test cases with 3000+ tubes have been calculated. 

Each tube has some properties like 3D position of start and end point, diameter, wall thickness, heat absorption, bend radius, orifice etc.

Setting up a data file for several 1000 tubes is quite time consuming. 

Therefore a kind of pre-processor is used. The name of the program is "dxf2wsc.exe". The whole input (as well as output) is based on geometry data in .dxf-format. Any program that can generate 3D-line dxf data and have the ability of a big number of layers with long layer names can be used. Unfortunately, most of open source design or drawing programs fail in one of the preconditions (AFAIK). 

The geometry has to be set up in the design program and each tube has to be on one layer with the additional data. Similar tubes can use the same layer. 

"dxf2wsc" reads the dxf-file and writes the data file for the circulation program. If there are problems with the dxf-file those problems will be given as error messages. If there is a problem with the geometry, like tube ends that are not connected to other tubes, a dxf-file will be produced. It can be overlayed with the original drawing and indicates the position of the problem. 

Finally, the general data like drum pressure, calculation method, tube roughness, convergence limit, maximum iterations etc. have to be given. 

In the second step the calculation program "wsc.exe" can be started. Here again errors from geometry like no path between downcomer and heated tube or no path between heated tube and drum can be found give the indication of the place as dxf-file. The calculation stops and the geometry should be corrected. If everything ok, the sum of flow difference in each iteration step appear on the screen. Those numbers should go down. 

After a successful calculation dxf-files for the different results are written. The most important one should be *projectname*_result_safety. If every line (tube) green -> everything ok. If the lines are red at least one safety criterion is violated and the circulation should be improved. If there are horizontal tubes also *projectname*_result_pattern should be checked for red lines. Certainly, only horizontal and heated tubes should be regarded.

This program is intended to calculate a stable operating condition. Load or pressure changes over a short time can also cause problems with the water circulation. Those cases are not covered by this program. 

### Bidrum boilers
The program is intended to calculate the water circulation in single drum boilers with unheated downcomers. However, it is even possible to calculate bidrum boilers, i.e. the downcomers are heated. In this case the tubes that should act as downcomer should be connected directly to steam drum and should have a given inlet enthalpy below saturated water enthalpy. Typically, in bidrum boilers with heated downcomers the feedwater pipe is adjacent to the downcomer inlet. The mixture enthalpy between drum water and feed water should assure that there is no evaporation in the downcomers. 

This should be regarded as experimental as there is no check if the inlet enthalpies are realistic. 
