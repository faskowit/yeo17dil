%% workflow to upsample, then dialte, the yeo 17 split annotation

% ./data/ directory with links to my local fsaverage and fsaverage5
% also ./data/yeoData_fsaverage5 with the data from the CBIG lab parcs

% first, using FreeSurfer functs, upsample the fsaverage5 annotations

% https://mail.nmr.mgh.harvard.edu/pipermail//freesurfer/2015-December/042811.html
% mri_surf2surf --srcsubject fsaverage5 --trgsubject fsaverage --hemi rh --sval-annot yeoData_fsaverage5/rh.Yeo2011_17Networks_N1000.split_components.annot --tval rh.yeo17splitUpSmp.annot
% mri_surf2surf --srcsubject fsaverage5 --trgsubject fsaverage --hemi lh --sval-annot yeoData_fsaverage5/lh.Yeo2011_17Networks_N1000.split_components.annot --tval lh.yeo17splitUpSmp.annot

% moved these outputs to folder:
% ./data/yeoData_fsaverage/

%% now time for dilation
% load data

addpath(genpath('~/JOSHSTUFF/projects/dilateParcellation/')) 
addpath(genpath('./data/'))

inf_surf_lh = './data/fsaverage/surf/lh.inflated' ;
inf_surf_rh = './data/fsaverage/surf/rh.inflated' ;

annot_lh = './data/yeoData_fsaverage/lh.yeo17splitUpSmp.annot' ;
annot_rh = './data/yeoData_fsaverage/rh.yeo17splitUpSmp.annot' ;

[~,annotLabs_lh] = read_annotation(annot_lh) ;
[~,annotLabs_rh,annotTable] = read_annotation(annot_rh) ;

[verts_lh,faces_lh] = read_surf(inf_surf_lh);
[verts_rh,faces_rh] = read_surf(inf_surf_rh);

%% neighbors

% left
s_lh = struct() ;
s_lh.nverts = size(verts_lh,1) ;
s_lh.nfaces = size(faces_lh,1) ;
s_lh.faces = faces_lh + 1 ;
s_lh.coords = verts_lh ;
n_lh = fs_find_neighbors(s_lh) ;

% right
s_rh = struct() ;
s_rh.nverts = size(verts_rh,1) ;
s_rh.nfaces = size(faces_rh,1) ;
s_rh.faces = faces_rh + 1 ;
s_rh.coords = verts_rh ;
n_rh = fs_find_neighbors(s_rh) ;

%% get vertex-wise data

w_lh = ones(length(annotLabs_lh),1);
w_rh = ones(length(annotLabs_rh),1);

for idx = 1:size(annotTable.table,1)
    disp(idx)
    w_lh(annotLabs_lh == annotTable.table(idx,5)) = idx;
    w_rh(annotLabs_rh == annotTable.table(idx,5)) = idx;
end

%% quick viz
% figure
% colormap(annotTable.table(:,1:3) ./ 255)
% 
% tmpStruct = s_rh ;
% tmpWei = w_rh ;
% 
% viewww = trisurf(tmpStruct.faces,...
%     tmpStruct.coords(:,1),...
%     tmpStruct.coords(:,2),...
%     tmpStruct.coords(:,3),...
%     tmpWei);
% set(viewww,'EdgeColor','none');
% axis equal; axis off
% view(90,0)
% camlight headlight; material dull; lighting gouraud
% viewww.CDataMapping = 'direct' ;

%% need medial wall too
setenv('SUBJECTS_DIR',strcat(pwd,'/data/'))
sname = 'fsaverage' ;

lname = 'lh.Medial_wall' ;
medial_wall_lh = read_label(sname,lname) ;
lname = 'rh.Medial_wall' ;
medial_wall_rh = read_label(sname,lname) ;

%% do the dilation

% left
w_dil_lh = dil_surf_parc(w_lh,n_lh.nbrs,1,2) ;
fin = 0 ;
idx = 1 ;
while ~fin
    disp(idx)
    idx = idx + 1;
    [w_dil_lh,fin] = dil_surf_parc(w_dil_lh,n_lh.nbrs,1,2) ;
end
% write over medial wall
w_dil_lh(medial_wall_lh(:,1)+1) = 1 ;

% right
w_dil_rh = dil_surf_parc(w_rh,n_rh.nbrs,1,2) ;
fin = 0 ;
idx = 1 ;
while ~fin
    disp(idx)
    idx = idx + 1;
    [w_dil_rh,fin] = dil_surf_parc(w_dil_rh,n_rh.nbrs,1,2) ;
end
% write over medial wall
w_dil_rh(medial_wall_rh(:,1)+1) = 1 ;

%% quick viz
% figure
% colormap(annotTable.table(:,1:3) ./ 255)
% 
% tmpStruct = s_rh ;
% tmpWei = w_dil_rh ;
% 
% viewww = trisurf(tmpStruct.faces,...
%     tmpStruct.coords(:,1),...
%     tmpStruct.coords(:,2),...
%     tmpStruct.coords(:,3),...
%     tmpWei);
% set(viewww,'EdgeColor','none');
% axis equal; axis off
% view(-90,0)
% camlight headlight; material dull; lighting gouraud
% viewww.CDataMapping = 'direct' ;

%% write annotation

%function write_annotation(filename, vertices, label, ct)
% vertices expected to be simply from 0 to number of vertices - 1;
% label is the vector of annotation
%
% ct is a struct
% ct.numEntries = number of Entries
% ct.orig_tab = name of original ct
% ct.struct_names = list of structure names (e.g. central sulcus and so on)
% ct.table = n x 5 matrix. 1st column is r, 2nd column is g, 3rd column
% is b, 4th column is flag, 5th column is resultant integer values
% calculated from r + g*2^8 + b*2^16 + flag*2^24. flag expected to be all 0

new_annotLabs_lh = ones(length(annotLabs_lh),1);
for idx = 1:size(annotTable.table,1)
    disp(idx)
    new_annotLabs_lh(w_dil_lh == idx) = annotTable.table(idx,5);
end

new_annotLabs_rh = ones(length(annotLabs_rh),1);
for idx = 1:size(annotTable.table,1)
    disp(idx)
    new_annotLabs_rh(w_dil_rh == idx) = annotTable.table(idx,5);
end

dil_annot_name = strcat(pwd,'/lh.yeo17dil.annot') ;
write_annotation(dil_annot_name,0:(length(new_annotLabs_lh)-1),new_annotLabs_lh,annotTable)

dil_annot_name = strcat(pwd,'/rh.yeo17dil.annot') ;
write_annotation(dil_annot_name,0:(length(new_annotLabs_rh)-1),new_annotLabs_rh,annotTable)
