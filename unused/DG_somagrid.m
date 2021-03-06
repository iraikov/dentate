%% Create Grid of Packed Somata and Select Those That Lie Within
%% the layer volume

function [Soma_Points] = DGnetwork_somagrid(SomaDistH,SomaDistV,layer_eq,layer_min,layer_mid,layer_max)

% Define parameters
Kepler      =   pi / (3 * sqrt (2));

% Create layer surface
[x_g1,y_g1,z_g1]    = feval(layer_eq,layer_min);
[x_g2,y_g2,z_g2]    = feval(layer_eq,layer_mid);
[x_g3,y_g3,z_g3]    = feval(layer_eq,layer_max);

X_g                 = [x_g1;x_g2;x_g3];
Y_g                 = [y_g1;y_g2;y_g3];
Z_g                 = [z_g1;z_g2;z_g3];
[~,S_g]             = alphavol([X_g(:),Y_g(:),Z_g(:)],120);

% Define limits for somata grid
xmin    =   -4000;
xmax    =   4000;
ymin    =   -500;
ymax    =   5000;
zmin    =   -900;
zmax    =   1800;

% Round limits to nearest multiple of soma size
xmax2    = round((xmax-xmin)/SomaDistH)*SomaDistH + xmin;
ymax2    = round((ymax-ymin)/SomaDistH)*SomaDistH + ymin;
zmax2    = round((zmax-zmin)/(Kepler*SomaDistV))*(Kepler*SomaDistV) + zmin;
zlayers  = ((zmax2-zmin)/(Kepler*SomaDistV));

% Create square packed somata within limits
xlin        = linspace(xmin,xmax2,(xmax2-xmin)/SomaDistH);
ylin        = linspace(ymin,ymax2,(ymax2-ymin)/SomaDistH);
zlin        = linspace(zmin,zmax2,zlayers);
[xx,yy,zz]  = meshgrid(xlin,ylin,zlin);
M           = [reshape(xx,[],1),reshape(yy,[],1),reshape(zz,[],1)];

% Shift every other z-layer to create hexagonal packing
xycells     = length(xlin)*length(ylin);
SomataGrid  = M(:,:);
for counter = 2:2:zlayers
    startpt                     = xycells*(counter-1) + 1;
    endpt                       = xycells*counter;
    SomataGrid(startpt:endpt,:) = [M(startpt:endpt,1)-SomaDistH,M(startpt:endpt,2),M(startpt:endpt,3)];
end

% Test whether somata lie within GCL and choose only those with centers inside GCL volume
IN_somata   = inpolyhedron(S_g.bnd,[X_g,Y_g,Z_g],SomataGrid,'FLIPNORMALS','True');
Soma_Points = SomataGrid(IN_somata,:);
