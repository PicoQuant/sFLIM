clc;
tic
[fname, pname] = uigetfile;
% fname = '14A_862Z_01_FRET_Slide_04_01.ptu';

name = [pname fname];
head = PTU_Read_Head(fname);

% [head, im_tcspc, im_chan, im_line, im_col, tag] = PTU_ScanRead(fname, 1, head);
[head, im_tcspc, im_chan, im_line, im_col, tag] = PTU_ScanRead(name, 1);

time_end  = toc/60