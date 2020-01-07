% available scenarios
% 4: 20 Devices with new assumptions 16/12/2019

function caseStudyData=callDatabase(scenario)

caseStudyData=[];

switch scenario
%     case 1
%          file_name='..\case_studies\Energies_paper.mat';
%                 load(file_name)
%     case 2
%          file_name='..\case_studies\Seven_paper.mat';
%                 load(file_name)
%                 
%     case 3
%            file_name='..\case_studies\Energies_paper20D.mat';
%                 load(file_name)
                
    case 4
          file_name='..\case_studies\Energies_paper20D_jan.mat';
                load(file_name)
        
end

end
