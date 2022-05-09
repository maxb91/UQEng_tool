function M_avg_strain = TensionChordModel(X)

% Fixed Parameters
D = 16; % Reinforcing bar diameter in mm
Es = 205000; % Elastic modulus of the steel in MPa
esu = 0.08; % Ultimate strain of the steel
fsu = 600; % Ultimate stress of the steel in MPa
fsy = 550; % Yield stress of the steel in MPa

% Material parameters for Ramberg-Osgood steel material law
ka = 0.002;
kb = 0.002;
alpha = log((esu-fsu/Es)/ka)/log(fsu/fsy);
kc =fsy/(kb^(1/alpha));

sigsr = X(:,1);
srm = X(:,2);
taub1 = X(:,3);

% Model - Alvarez 1998
% "Tension Chord Model" - calculate average strain in the reinforcing bar
M_avg_strain = sigsr/Es-taub1.*srm/(Es*D)+D./(2*taub1.*srm)*1./(alpha+1)./kc.^alpha.*(sigsr.^(alpha+1)-(sigsr-2*taub1.*srm/D).^(alpha+1));