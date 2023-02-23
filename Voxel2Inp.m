%% Voxel to INP
clear;clc;fclose all;close all;
% LOAD PROCESSED STRUTS
load('FCCZ_0.3_5x5x5_2x2x2_strutsProcessed.mat');

latticeTop='FCCZ_25_11_21';
% MAKE A DIRECTORY FOR THE SPECIFIC STRUT
mkdir(latticeTop)

%FCC CT RESOLUTION 0.010000600000000mm %FCCZ CT RESOLUTION 0.010000600000000
ctRes=0.010000600000000; % mm
for i=1:3:length(strutsProcessed)%length(input_file.need_to_rerun)randperm(length(strutsProcessed),150)
    %% GET NODES AND ELEMENTS
    %CREATE SPECIFIC STRUT FOLDER
    strut_id=i;
    fldrName=sprintf('%s\\%s_Strut_%d',latticeTop,latticeTop,strut_id);
    mkdir(fldrName);
    
    filenameINP=sprintf('%s\\%s_Strut_%d.inp',fldrName,latticeTop,strut_id);
    fileID=fopen(filenameINP,'w');
    binaryStack=strutsProcessed{i}; % strut number 1 of 1200 in strutsprocessed
%     volshow(strutsProcessed{i},'Renderer','VolumeRendering');


    [FEnodes,elems]=generateMesh(binaryStack);
    %% Get Nodal Centres
    FEnodes=[FEnodes(:,1) (FEnodes(:,2:end)-mean(FEnodes(:,2:end),1))*ctRes];
    
    % Centroid
    zsorted=unique(FEnodes(:,4));
    zupper=zsorted(end-1:end);zlower=zsorted(1:2);
    % Upper IDs
    upperID=FEnodes(:,4)==zupper(1) | FEnodes(:,4)==zupper(2);
    % Lower IDs
    lowerID=FEnodes(:,4)==zlower(1) | FEnodes(:,4)==zlower(2);
    
  
    %% WRITE STRUT INP
    printHeaders(fileID);
    fprintf(fileID,'*Node\n');
    fprintf(fileID,'      %d,        %0.3f,        %0.3f,          %0.3f\n',...
        [FEnodes(:,1) FEnodes(:,2) FEnodes(:,3) FEnodes(:,4)]');
    fprintf(fileID,'*Element, type=C3D8R\n');
    fprintf(fileID,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',elems');
    fprintf(fileID,'\n*Elset, elset=Set-1, generate\n');
    fprintf(fileID,'\t1,  %d,       1\n',max(elems(:,1)));
    fprintf(fileID,'** Section: Section-1\n');
    fprintf(fileID,'*Solid Section, elset=Set-1, material=TI6AL4V\n');
    fprintf(fileID,'\n,\n');
    fprintf(fileID,'*End Part\n**\n**\n');
    fprintf(fileID,'** ASSEMBLY\n**\n');
    fprintf(fileID,'*Assembly, name=Assembly\n**\n');
    fprintf(fileID,'*Instance, name=PART-1-1, part=PART-1\n*End Instance\n**\n');

    % Upper And Lower Reference Nodesets
    % Calculate Centroids
    upperRefXYZ=[mean(FEnodes(upperID,2)) mean(FEnodes(upperID,3)) zupper(2)+(5*ctRes)];
    lowerRefXYZ=[mean(FEnodes(lowerID,2)) mean(FEnodes(lowerID,3)) zlower(1)-(5*ctRes)];
    fprintf(fileID,'*Node\n');
    fprintf(fileID,'1,           %0.3f,           %0.3f,         %0.3f\n',upperRefXYZ(1),upperRefXYZ(2),upperRefXYZ(3));
    fprintf(fileID,'2,           %0.3f,           %0.3f,         %0.3f\n',lowerRefXYZ(1),lowerRefXYZ(2),lowerRefXYZ(3));
    fprintf(fileID,'*Nset, nset=UpperNodeSet\n');
    fprintf(fileID,' 1,\n');
    fprintf(fileID,'*Nset, nset=LowerNodeSet\n');
    fprintf(fileID,' 2,\n');
    % Reference Sets
    nodeSetUpper=FEnodes(upperID,1);
    id_upper=floor(length(nodeSetUpper)./16);
    nodeSetLower=FEnodes(lowerID,1);
    id_lower=floor(length(nodeSetLower)./16);
    
    
    % Upper
    fprintf(fileID,'*Nset, nset=UpperNodeSurfaceSet, instance=PART-1-1\n');
    fprintf(fileID,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',nodeSetUpper(1:id_upper*16));
    fprintf(fileID,'%d, ',nodeSetUpper(id_upper*16+1:end-1));
    fprintf(fileID,'%d\n',nodeSetUpper(end));
    % Lower
    fprintf(fileID,'*Nset, nset=LowerNodeSurfaceSet, instance=PART-1-1\n');
    fprintf(fileID,'%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',nodeSetLower(1:id_lower*16));
    fprintf(fileID,'%d, ',nodeSetLower(id_lower*16+1:end-1));
    fprintf(fileID,'%d\n',nodeSetLower(end));
    
    %Surface Defintion
    fprintf(fileID,'*Surface, type=NODE, name=UpperSurf\n');
    fprintf(fileID,'UpperNodeSurfaceSet, 1.\n');
    fprintf(fileID,'*Surface, type=NODE, name=LowerSurf\n');
    fprintf(fileID,'LowerNodeSurfaceSet, 1.\n');
    
    % Coupling Constraint
    fprintf(fileID,'** Constraint: Constraint-1\n');
    fprintf(fileID,'*Coupling, constraint name=Constraint-1, ref node=UpperNodeSet, surface=UpperSurf\n');
    fprintf(fileID,'*Kinematic\n');
    fprintf(fileID,'** Constraint: Constraint-2\n');
    fprintf(fileID,'*Coupling, constraint name=Constraint-2, ref node=LowerNodeSet, surface=LowerSurf\n');
    fprintf(fileID,'*Kinematic\n');
    fprintf(fileID,'*End Assembly\n');

    
    
    %% MATERIAL PROPERTIES
    fprintf(fileID,'**\n** MATERIALS\n**\n');
    fprintf(fileID,'*Material, name=TI6AL4V\n');
    fprintf(fileID,'*Density\n');
    fprintf(fileID,'4.43e-09,\n');
    fprintf(fileID,'*Elastic\n');
    fprintf(fileID,'113800., 0.342\n');
    
    %% STEP DEF
    fprintf(fileID,'** ----------------------------------------------------------------\n**\n');
    fprintf(fileID,'** STEP: Step-1\n**\n');
    fprintf(fileID,'*Step, name=Step-1, nlgeom=NO, perturbation\n');
    fprintf(fileID,'*Buckle\n');
    fprintf(fileID,'5, , 10, 30\n**\n');
    
    %% Boundary Conditions
    fprintf(fileID,'** BOUNDARY CONDITIONS\n**\n');
    fprintf(fileID,'** Name: BC-1 Type: Symmetry/Antisymmetry/Encastre\n');
    fprintf(fileID,'*Boundary, op=NEW, load case=1\n');
    fprintf(fileID,'LowerNodeSet, ENCASTRE\n');
    fprintf(fileID,'*Boundary, op=NEW, load case=2\n');
    fprintf(fileID,'LowerNodeSet, ENCASTRE\n**\n');
    % UNIT LOAD
    fprintf(fileID,'** LOADS\n**\n');
    fprintf(fileID,'** Name: Load-1   Type: Concentrated force\n');
    fprintf(fileID,'*Cload\n');
    fprintf(fileID,'UpperNodeSet, 3, -1.\n');
    
    fprintf(fileID,'**\n** OUTPUT REQUESTS\n**\n');
    fprintf(fileID,'*Restart, write, frequency=0\n');
    fprintf(fileID,'**\n** FIELD OUTPUT: F-Output-1\n**\n');
    fprintf(fileID,'*Output, field, variable=PRESELECT\n');
    fprintf(fileID,'*End Step');
    
    fclose(fileID);
    
    %% RUN THE JOB
    folderPathFinal=fullfile(pwd,fldrName);
    % DEPENDING ON WINDOWS VERSION & or | must be used to sepreate the CD
    % line from the input line
    cmd=sprintf('cd "%s" & abaqus interactive job=%s_Strut_%d cpus=16',folderPathFinal,latticeTop,strut_id);
    [status, result] = dos(cmd);

end



%% Functions
function printHeaders(fileID)
fprintf(fileID,'*Heading\n');
fprintf(fileID,'** Job name: Job-1 Model name: Model-1\n');
fprintf(fileID,'** Generated by: Abaqus/CAE 6.14-1\n');
fprintf(fileID,'*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');
fprintf(fileID,'**\n');
fprintf(fileID,'** PARTS\n');
fprintf(fileID,'**\n');
fprintf(fileID,'*Part, name=Part-1\n');
end
function [FEnodes,elems]=generateMesh(binaryStack)
testImg=binaryStack(:,:,1);
[~,pPerCol]=size(testImg);
gridPbaseLeft=@(i,j,total,iter) ((pPerCol+1)*i-pPerCol+j-1)+(iter-1)*total;
%% Test
xrange=0.5:1:size(testImg,1)+0.5;
yrange=0.5:1:size(testImg,2)+0.5;
zrange=1:size(binaryStack,3)+1;
imTotalPts=length(xrange)*length(yrange);
%%
[X,Y,Z]=meshgrid(xrange,yrange,zrange);
nodes=[X(:),Y(:),Z(:)];
numElemsPerLayer=arrayfun(@(x) sum(binaryStack(:,:,x)==1,'all'),1:size(binaryStack,3));
totalNumOfElements=sum(numElemsPerLayer);
elemArray=cell(size(binaryStack,3),1);
for i=1:size(binaryStack,3)
    %%Find binary ims
    [x,y]=find(binaryStack(:,:,i)==1);
    vals1=arrayfun(@(x,y) [(gridPbaseLeft(x,y,imTotalPts,i)+pPerCol+1),(gridPbaseLeft(x,y,imTotalPts,i)+pPerCol+2),...
        gridPbaseLeft(x,y,imTotalPts,i)+1 gridPbaseLeft(x,y,imTotalPts,i)],x,y,'UniformOutput',0);
    vals2=cellfun(@(x) x+imTotalPts,vals1,'UniformOutput',0);
    elemArray{i}=cell2mat([vals1 vals2]);
end
%%
elemArray=cell2mat(elemArray);
[C,~,~] = unique(elemArray);
nodesUpdated=nodes(C,:);
FEnodes=[C nodesUpdated];
elems=[(1:totalNumOfElements)' elemArray];

end

