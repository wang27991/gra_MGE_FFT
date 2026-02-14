# A-Program-of-the-Forward-Modeling-of-MGE-and-the-Inversion-Using-FFT
This is a software provides a complete solution for gravity data inversion using the Fast Forward Modeling of MGE method and the Inversion method Using FFT. The suite includes programs for fast gravity kernel matrix computation and efficient 3D density inversion with FFT.
How do I get set up? Copy the repository folder (*** Downloads *** from the menu on the left of this page) in Windows 64bit. In alternative use any other machine with matlab 2019b or later installed.

Usage

Choose one of the following:
% Open and modify A_the MGE.m bbb = 'your_data_identifier'; % Must match your data files nz = 20; % Vertical grid cells zmin = -3000; % Model bottom depth (meters)

% Run preprocessing run('A_the MGE.m');

This generates GgF{identifier}.mat containing the Fourier-transformed kernel matrix.

% Open and modify B_main_inv_fft.m bbb = 'your_data_identifier'; % Must match preprocessing zmin = -10000; % Inversion model bottom depth Itermax = 100; % Maximum iterations
% Run inversion run('B_main_inv_fft.m');

OUTPUT FILES :

.mat> containing the Fourier-transformed kernel matrix. <gravinv30_20_10_inv>.txt --> 3D inversion result

Who do I talk to? Guoqing Ma, Jilin University *** maguoqing@jlu.edu.cn


tips:All the test files have been provided in the repository. Download all the files into one package and run them in order as per the instructions in the Readme.
