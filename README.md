# Boiler Water Circulation ![Here should be the logo](WSCLogo.png)

The program calculates the flow in the network of tubes that represent a boiler. It calculates the densities (as driving force) and resistances. The resulting flow will be compared to some criteria for safe flow. 

The calculation is done on a network of tubes limited only by the storage capacity of the computer. Test cases with 3000+ tubes have been calculated. 

Each tube has some properties like 3D position of start and end point, diameter, wall thickness, heat absorption, bend radius, orifice etc. Those data are contained in a .dxf-file, that can be generated with a CAD program. Any program that can generate 3D-line .dxf data and have the ability of a big number of layers with long layer names can be used.

A kind of pre-processor is used. The name of the program is "dxf2wsc.exe". It reads the dxf-file and writes the data file for the circulation program. If there are problems with the dxf-file those problems will be given as error messages and a .dxf-file with lines pointing to the line (tube) or point where this error occured. 

In the second step the calculation program "wsc.exe" can be started. After a successful calculation dxf-files for the different results like flow or safety criteria are written. 

Transient operation, like load or pressure changes, are not covered by this program.  

### Getting Started

Download the files "dxf2wsc.exe" and "wsc.exe" if you want to run this program under Windows. Copy it to one folder. A folder dedicated to this calculation would be a good idea. Set up the underlying file system as described in [mainpage] (/wsc2dxf/doxygen/mainpage.md). 

If you want to run the programs under a different operation system, you can compile it from source. Both programs are "console" programs. There is no specific windowing system involved. 

### Prerequisites
To clone the repositories you should have git (https://git-scm.com/) installed or have git-add-on to your preferred IDE activated. 

The preprocessor "dxf2wsc" only needs a recent C++ compiler to be compiled and linked. Both programs use STL functionality (C++ 11 and following), mainly vectors and lists.

The main program "wsc" also needs a recent C++ compiler. MinGW64(https://www.mingw-w64.org/) , clang(https://clang.llvm.org/) and VisualC(https://visualstudio.microsoft.com/) work well with this program. *eigen* (https://eigen.tuxfamily.org/index.php?title=Main_Page) is used for solving the linear equation system.

*eigen* comes as a set of header files. Download *eigen* to a folder of your choice. "additional header files" should point to this folder in your make-file. The program was testet with version 3.3 of *eigen* (a compilation with *eigen3.4* leads to a lot of errors within *eigen*)

The programs can be compiled using a make-file only or using an integrated development environment (IDE). During the development Netbeans (https://netbeans.apache.org/), Code::Blocks(https://www.codeblocks.org/), KDevelop(https://www.kdevelop.org/) and Visual Studio Community (https://visualstudio.microsoft.com/) were used. In the moment VSCode(https://code.visualstudio.com/) is quite popular, never tried it. Running VSCode inside the "source" folder should work but never tried. 

### Installing

As said above there are different ways to compile the program. I would use an IDE and make a new empty project in a folder of my choice. Close the IDE. Run GIT, clone the repository and checkout. Alternatively you can clone directly from Github to this folder and checkout. Start the IDE again with the empty project. Add all files in "source" folder to project. Adjust the build configuration to your compiler and linker. A 64 bit configuration with minimum C++ 11 standard is recommended. Add the folder that contains *eigen* to your header search path. 

Build.

Copy the resulting "wsc.exe" to your working directory (see above)

"dxf2wsc" doesn't rely on *eigen*. Therefore it's not needed in the search path.  

### Running the tests

Copy the testfile (testfile.dxf and testfile.readme) to (your_folder_name\WSCDATA\testfile\testfile.dxf). Open a console, f.i. cmd or PowerShell. Run "dxf2wsc.exe" and give the project name "testfile". Additional data like pressure etc. can be found in "testfile.readme". 

After completition run "wsc.exe" and give project name again. "wsc" should terminate without errors. Now you can look at the results "your_folder_name\WSCDATA\testfile\testfile_result_***.dxf". Most important should be "testfile_result_safety.dxf". Open it in your CAD program. You will see the tube (line) structure as in original setup. If the lines are in green -> everything fine. If a line is red we have to expect some problem. Hovering with the mouse over the line will show the layer name, which contains detail information like safety factor.


### Built With

* [*eigen*] (https://eigen.tuxfamily.org/index.php?title=Main_Page) - solvers for linear equation systems 
* [git] (https://git-scm.com/) - for cloning repository
* [doxygen] (https://www.doxygen.nl/index.html) - for generation of documentation

### Contributing

Please read [CONTRIBUTING.md](https://github.com/WizenedEngineering/BoilerWaterCirculation/code_of_conduct.md) for details on our code of conduct, and the process for submitting pull requests to us.

### Versioning

For the versions available, see the [tags on this repository](https://github.com/WizenedEngineering/BoilerWaterCirculation/tags). 

### Authors

* **Rainer Jordan** - *Initial work* - [WizenedEngineering](https://github.com/WizenedEngineering)

See also the list of [contributors](https://github.com/WizenedEngineering/BoilerWaterCirculation/contributors) who participated in this project.

### License

This project is licensed under the European Union Public Licence (EUPL), Version 1.2 or later - see the [LICENSE.md](LICENSE.md) file for details

### Acknowledgments

* Thanks to Michael Beyer and Stefan Kohn for the valuable discussion

