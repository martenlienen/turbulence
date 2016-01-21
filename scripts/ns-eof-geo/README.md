Programm for creating an arbitrary geometry

Steps:
1) Open 'mainbatch.m' and set reference to the config file (it contains all informations of geometry and mesh)
2) Choose and specify information to the arbitrary geometry
3) Run 'mainbatch.m'

Bug:
Before pasing the config file, you need to modify slightly the following line:
	<mesh t="uniform">uniform</mesh>
This is necessary because MATLAB seems to have to access the XML-DOM correctly and cannot read the node value.
