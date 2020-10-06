This folder contains code to run the model from 

Chariker L, Shapley R, Young LS (2016) Orientation selectivity from very sparse LGN inputs in a comprehensive model of macaque V1 cortex. J Neurosci 36:12368 –12384.

The model is coded in Java (primarily in the file Network.java), and it can be loaded and run from within Matlab.  This readme contains basic instructions for compiling the Java code using the Netbeans IDE, and for running the accompanying Matlab scripts to show model output.

COMPILING JAVA CODE

1. Download and install the Java Development Kit with NetBeans IDE, which can be found at

https://www.oracle.com/technetwork/java/javase/downloads/jdk-netbeans-jsp-3413139-esa.html

for Linux, Windows, or Mac OS X.

2. Run the NetBeans IDE.  In the File drop down menu click "New Project...".  Select "Java Application" from the project types and click "Next".  In the Project Name field, type "Sim" and click "Finish".  A new project should be created with the name "Sim" and a single package named "sim".  In some versions of Netbeans, the package will be named "main".  In that case, you can right-click on the package, click on "Refactor->Rename", and rename the package "sim".

3. Next, because many versions of Matlab are not compatible with the latest version of Java (version 8 or later), you will need to compile the project with version 6 or 7.  To do this in NetBeans, right click on the project name, and click "Properties".  From there, in the "Source/Binary Format" field, select "JDK 7", and click "OK".

4. Next, two third-party packages are to be added to the project.  Under the "Sim" project, right-click on the "Libraries" folder and select "Add JAR/folder".  Navigate to this zip file and select the .jar files "colt-1.2.0" and "jmatio", select them, and click "Open".  This should add the packages to the project.

5. Add all .java files from this folder to the "sim" package.  The files can be dragged and dropped from this folder into the package in NetBeans.

6. Finally, click on "Clean and Build Project".  If the code compiles correctly, the output window should print the line "BUILD SUCCESSFUL".

RUNNING ACCOMPANYING MATLAB SCRIPTS

The Network class, defined by "Network.java", contains all model parameters as well as methods for building and simulating the model.  The loadjava.m Matlab script loads an instance of the Java class "Network" into memory with object name "a".  Model parameters can then be assigned by accessing the appropriate fields of "a".  loadjava.m also runs methods for allocating memory and building the network connectivity.  

In order to run loadjava.m, three lines of the script must be edited to include paths to the compiled Java code as well as two .mat.  The comments in the script indicate where the paths must be placed.  The .mat files are in this folder, and the compiled Java code is located in the output folder (likely named "dist") of the NetBeans project directory.  Once the paths are input, the script can be run, and the model loaded into memory.

The script "test.m" can then be run, which initiates model parameters, runs the model with background drive, and outputs the E and I firing rates (in red and blue, respectively) over time in two patches (the vertical-preferring patch located at coordinates (500,750), and the horizontal-preferring patch located at (1000,750)).





