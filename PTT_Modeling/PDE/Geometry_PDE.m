%% -------------------------------PDE Solution------------------------------
clear;clc;
skin = [3;4;0;0.13;0.13;0;3;3;-3;-3];
tumor = [1;0;0;0.5];
rect2 = [3;4;-0.5;0;0;-0.5;0.5;0.5;-0.5;-0.5];
adipose = [3;4;0.13;0.3;0.3;0.13;3;3;-3;-3];
muscle = [3;4;0.3;2;2;0.3;3;3;-3;-3];

tumor = [tumor;zeros(length(skin) - length(tumor),1)];
gd = [skin,tumor,rect2,adipose,muscle];

ns = char('skin','tumor','rect2','adipose','muscle');
ns = ns';

sf = '(skin+adipose+muscle+tumor)-rect2';

[dl,bt] = decsg(gd,sf,ns);

pdegplot(dl,'EdgeLabels','on','FaceLabels','on')
xlim([-1,2.2])
axis equal