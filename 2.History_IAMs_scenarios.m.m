%% Recalculate actual_regions_emissions and compare models data with history

clear all

[~,~,CO2_ISO_total]=xlsread('F:\research\code_and_data\data_generation\ISO_ar6_mainmodel.csv');
[~,~,CO2_R5_total]=xlsread('F:\research\code_and_data\data_generation\R5_ar6_mainmodel.csv');
CO2_total_cha=[CO2_ISO_total;CO2_R5_total(2:end,:)];
clearvars -except CO2_total_cha

[~,~,CO2_hist_own]=xlsread('F:\research\code_and_data\data_generation\Comparison of emissions for 60 developing countries with IEA.xlsx'); 
[~,~,CO2_hist_IEA]=xlsread('F:\research\code_and_data\data_generation\institution_comparison.xlsx','IEA_MT'); 
CO2_hist_own_summedMatrix=[CO2_hist_own(:,1:3);CO2_hist_own(2:end,5:7)];
CO2_hist_IEA=transpose(CO2_hist_IEA);

% === Format historical own data to standard format with ISO code ===
CO2_h_final_own(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','own_ISO'); 
CO2_coun_unique = unique(CO2_hist_own_summedMatrix(:,1));
for i=1:length(CO2_coun_unique)
    CO2_1=CO2_hist_own_summedMatrix(find(ismember(CO2_hist_own_summedMatrix(:,1),CO2_coun_unique(i,:))==1),:);
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(:,3)'];
    if length(aa)==13
    CO2_h_final_own=[CO2_h_final_own;aa];
        
    end
end
CO2_h_final_own(1,1:2)={'Country','ISO'};
CO2_h_final_own=CO2_h_final_own(:,1:13);

% === Format IEA data to standard format with ISO code ===
CO2_h_final_IEA(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','IEA_ISO'); 
error_index = strcmp(ISO_name_match(:, 2), 'ActiveX VT_ERROR: ');


ISO_name_match = ISO_name_match(~error_index, :);

CO2_hist_IEA=CO2_hist_IEA(2:end,:);
for i=1:length(CO2_hist_IEA)
    CO2_1=[CO2_hist_IEA(i,1),CO2_hist_IEA(i,12:22)];
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(1,2:end)];
    if length(aa)==13
    CO2_h_final_IEA=[CO2_h_final_IEA;aa];
        
    end
end
CO2_h_final_IEA(1,1:2)={'Country','ISO'};


% ===Combine own compiled emissions with IEA emissions===
CO2_h_final=[];
CO2_h_final_IEA=CO2_h_final_IEA(2:end,:);
for i=1:size(CO2_h_final_IEA,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1),:);
    else
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),CO2_h_final_IEA(i,2))==1),:);
    
    end
    
    CO2_h_final=[CO2_h_final;A];
end
% === Aggregate regional emissions using ISO-to-region mapping ===
 [mapping_country_model, ~, country_match] = xlsread('F:\research\code_and_data\data_generation\IAM_All_Country_Match(test) R5.xlsx','AIM|CGE');  
  mapping_country_model(isnan(mapping_country_model))=0;
  country_match=country_match(2:end,:);
  country_match(3:end-2,5:end)=num2cell(mapping_country_model);
  country_match(end-1:end,5:end)=num2cell(zeros(2, size(mapping_country_model, 2)));
  country_match=[country_match(:,1:9),country_match(:,12)];
  mapping_country_model=cell2mat(country_match(3:end,5:end));
 % Aggregate historical emissions by region
  region_aggregate=[];
   for i=5:size(country_match,2) 
        name_region=country_match(:,i);
        name_V1=country_match(3:end,1:2);
        indivmodel_name_region=name_V1(find(cell2mat(name_region(3:end,1))==1),:);
     carbon_R=[];
             for h=1:size(indivmodel_name_region,1)
            match_idx = find(strcmp(CO2_h_final(:,2), indivmodel_name_region(h,2)));
             CO2_c=[CO2_h_final(match_idx,2),CO2_h_final(match_idx,3:end)];
             carbon_R=[carbon_R;CO2_c];
             end 
             carbon_R1=sum(cell2mat(carbon_R(:,2:end)),1);
             hu=[name_region(1,1),name_region(1,1),num2cell(carbon_R1)];
              region_aggregate=[region_aggregate;hu];
   end    
CO2_h_final=[CO2_h_final;region_aggregate];
clearvars -except CO2_total_cha CO2_h_final
% === Match modelled data with historical emissions ===
CO2_actual_final=[];
CO2_origin_final=[];
for i=2:size(CO2_total_cha,1)
    hu=CO2_h_final(find(ismember(CO2_h_final(:,2),CO2_total_cha(i,8))==1),3:13);
    if ~isempty(hu)
    hui=[CO2_total_cha(i,1:10),hu];
    CO2_origin_final=[CO2_origin_final;CO2_total_cha(i,:)];
    CO2_actual_final=[CO2_actual_final;hui];
    
    end
end
title=[CO2_total_cha(1,1:10),num2cell((2010:1:2020))];
CO2_actual_final=[title;CO2_actual_final];
CO2_origin_final=[CO2_total_cha(1,:);CO2_origin_final];
xlswrite('F:\research\code_and_data\data_generation\emissions_total.xlsx',CO2_actual_final,'actual');
xlswrite('F:\research\code_and_data\data_generation\emissions_total.xlsx',CO2_origin_final,'origin_models');
%% all models all scenarios_average annual growth rate
clear all

[~,~,CO2_ISO_total]=xlsread('F:\research\code_and_data\data_generation\ISO_ar6_mainmodel.csv');
[~,~,CO2_R5_total]=xlsread('F:\research\code_and_data\data_generation\R5_ar6_mainmodel.csv');
CO2_total_cha=[CO2_ISO_total;CO2_R5_total(2:end,:)];
CO2_total=CO2_total_cha;
CO2_total=CO2_total(:,1:13);
[num_rows, num_cols] = size(CO2_total);  
aagr_vec = NaN(num_rows, 1);
cagr_vec = NaN(num_rows, 1);

% === Calculate 5-year average growth rates (2010–2015, 2015–2020) ===
for i = 2:num_rows
    row_data = CO2_total(i, 11:end);  
    vals = cellfun(@(x) str2double(string(x)), row_data);
    for j=1:2
     if ~isnan(vals(j)) && ~isnan(vals(j+1)) && vals(j) ~= 0
            aagr_cell{i,j} = (vals(j+1)/vals(j))^(1/5) - 1;
        else
            aagr_cell{i,j} = NaN;  
        end
    end
end
CO2_FINAL_1=[CO2_total(2:end,1:10),aagr_cell(2:end,:)];
CO2_FINAL=[CO2_total(1,1:end-1);CO2_FINAL_1];
clearvars -except CO2_FINAL

%% actual_average annual growth rate
clear all
[~,~,CO2_hist_own]=xlsread('F:\research\code_and_data\data_generation\Comparison of emissions for 60 developing countries with IEA.xlsx'); 
[~,~,CO2_hist_IEA]=xlsread('F:\research\code_and_data\data_generation\institution_comparison.xlsx','IEA_MT'); 
CO2_hist_own_summedMatrix=[CO2_hist_own(:,1:3);CO2_hist_own(2:end,5:7)];

CO2_hist_IEA=transpose(CO2_hist_IEA);

CO2_h_final_own(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','own_ISO'); 
CO2_coun_unique = unique(CO2_hist_own_summedMatrix(:,1));
for i=1:length(CO2_coun_unique)
    CO2_1=CO2_hist_own_summedMatrix(find(ismember(CO2_hist_own_summedMatrix(:,1),CO2_coun_unique(i,:))==1),:);
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(:,3)'];
    if length(aa)==13
    CO2_h_final_own=[CO2_h_final_own;aa];
        
    end
end
CO2_h_final_own(1,1:2)={'Country','ISO'};
CO2_h_final_own=CO2_h_final_own(:,1:13);

CO2_h_final_IEA(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','IEA_ISO'); 
error_index = strcmp(ISO_name_match(:, 2), 'ActiveX VT_ERROR: ');


ISO_name_match = ISO_name_match(~error_index, :);

CO2_hist_IEA=CO2_hist_IEA(2:end,:);
for i=1:length(CO2_hist_IEA)
    CO2_1=[CO2_hist_IEA(i,1),CO2_hist_IEA(i,12:22)];
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(1,2:end)];
    if length(aa)==13
    CO2_h_final_IEA=[CO2_h_final_IEA;aa];
        
    end
end
CO2_h_final_IEA(1,1:2)={'Country','ISO'};

% === Merge own compiled emissions with IEA emissions by ISO Code ===
CO2_h_final=[];
CO2_h_final_IEA=CO2_h_final_IEA(2:end,:);
for i=1:size(CO2_h_final_IEA,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1),:);
    else
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),CO2_h_final_IEA(i,2))==1),:);
    
    end
    
    CO2_h_final=[CO2_h_final;A];
end
% === Aggregate national emissions to regional level ===
 [mapping_country_model, ~, country_match] = xlsread('F:\research\code_and_data\data_generation\IAM_All_Country_Match(test) R5.xlsx','AIM|CGE');  
  mapping_country_model(isnan(mapping_country_model))=0;
  country_match=country_match(2:end,:);
  country_match(3:end-2,5:end)=num2cell(mapping_country_model);
  country_match(end-1:end,5:end)=num2cell(zeros(2, size(mapping_country_model, 2)));
  country_match=[country_match(:,1:9),country_match(:,12)];
  mapping_country_model=cell2mat(country_match(3:end,5:end));
  region_aggregate=[];
   for i=5:size(country_match,2) 
        name_region=country_match(:,i);
        name_V1=country_match(3:end,1:2);
        indivmodel_name_region=name_V1(find(cell2mat(name_region(3:end,1))==1),:);
     carbon_R=[];
             for h=1:size(indivmodel_name_region,1)
            match_idx = find(strcmp(CO2_h_final(:,2), indivmodel_name_region(h,2)));
             CO2_c=[CO2_h_final(match_idx,2),CO2_h_final(match_idx,3:end)];
             carbon_R=[carbon_R;CO2_c];
             end 
             carbon_R1=sum(cell2mat(carbon_R(:,2:end)),1);
             hu=[name_region(1,1),name_region(1,1),num2cell(carbon_R1)];
              region_aggregate=[region_aggregate;hu];
   end    
CO2_h_final=[CO2_h_final;region_aggregate];
clearvars -except CO2_h_final
% === Average Annual Growth Rates ===
average_annual_change_rate_1 = (cell2mat(CO2_h_final(:,8)) ./ cell2mat(CO2_h_final(:,3))).^(1/5) - 1;
average_annual_change_rate_2 = (cell2mat(CO2_h_final(:,13)) ./ cell2mat(CO2_h_final(:,8))).^(1/5) - 1;
average_annual_change_rate=[CO2_h_final(:,1:2),num2cell(average_annual_change_rate_1),num2cell(average_annual_change_rate_2)];

clearvars -except CO2_h_final average_annual_change_rate 

xlswrite('F:\research\code_and_data\data_generation\actual_change_rate.xlsx',average_annual_change_rate,'actual');
%%  downscaling_individual_countries_scenarios
clear all
cd 'F:\research\code_and_data\data_generation'
filename_1 = 'ISO_ar6_mainmodel.csv';  
CO2_ISO_total = readtable(filename_1, 'ReadVariableNames', true);
filename_2 = 'R5_ar6_mainmodel.csv';  
CO2_R5_total = readtable(filename_2, 'ReadVariableNames', true);
CO2_total_cha=[CO2_ISO_total;CO2_R5_total(2:end,:)];
% === Interpolation Module: Expand CO2 data from 5-year intervals to 1-year intervals ===
years = 2010:5:2100;
n_countries = height(CO2_total_cha);
result = array2table([years'], 'VariableNames', {'Year'});
CO2_total_cha_final=CO2_total_cha(:,1:10);
new_years = 2010:2100;
new_year_names = strcat("Y", string(new_years));
CO2_total_cha_final(:,11:101) = array2table(nan(height(CO2_total_cha), numel(new_years)));  
CO2_total_cha_final.Properties.VariableNames(11:101) = new_year_names;

for i = 1:n_countries
   
    
 
    row_data = table2array(CO2_total_cha(i, 11:end)); 
    year_idx = years(~isnan(row_data)); 
    values = row_data(~isnan(row_data));
    
   
    
    
    % === interp ===
    interp_years = 2010:1:max(year_idx);
    interp_values = pchip(year_idx, values, interp_years);
    CO2_total_cha_final(i,11:end)=array2table(interp_values);
   
    
end
% === Convert Table to Cell Format ===
% Extract variable names as a cell row
varNames = CO2_total_cha_final.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "Y");
dataCell = table2cell(CO2_total_cha_final);
CO2_total_cha_final_cell = [varNames; dataCell];
CO2_total_cha=CO2_total_cha_final_cell;


varNames = CO2_ISO_total.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "x");
dataCell = table2cell(CO2_ISO_total);
CO2_ISO_total = [varNames; dataCell];


varNames = CO2_R5_total.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "x");
dataCell = table2cell(CO2_R5_total);
CO2_R5_total = [varNames; dataCell];
clearvars -except CO2_total_cha CO2_ISO_total CO2_R5_total


models_name = CO2_R5_total(2:end,1);
models_name_1 = CO2_total_cha(:,1);

unique_models_name = unique(models_name);

% === Group Data by Model ===
% For each unique model, extract all matching rows from CO2_total_cha and store in cell array
model_data = cell(size(unique_models_name));
model_data_name=[];
for i = 1:numel(unique_models_name)
    model = unique_models_name(i);
    model_data{i}=CO2_total_cha(find(ismember(models_name_1,model)==1),:);
    model_data_name=[model_data_name;model];
end

[~,~,CO2_hist_own]=xlsread('F:\research\code_and_data\data_generation\Comparison of emissions for 60 developing countries with IEA.xlsx'); 
[~,~,CO2_hist_IEA]=xlsread('F:\research\code_and_data\data_generation\institution_comparison.xlsx','IEA_MT'); 
CO2_hist_own_summedMatrix=[CO2_hist_own(:,1:3);CO2_hist_own(2:end,5:7)];
CO2_hist_IEA=transpose(CO2_hist_IEA);

% === Format historical own data to standard format with ISO code ===
CO2_h_final_own(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','own_ISO'); 
CO2_coun_unique = unique(CO2_hist_own_summedMatrix(:,1));
for i=1:length(CO2_coun_unique)
    CO2_1=CO2_hist_own_summedMatrix(find(ismember(CO2_hist_own_summedMatrix(:,1),CO2_coun_unique(i,:))==1),:);
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(:,3)'];
    if length(aa)==13
    CO2_h_final_own=[CO2_h_final_own;aa];
        
    end
end
CO2_h_final_own(1,1:2)={'Country','ISO'};
CO2_h_final_own=CO2_h_final_own(:,1:13);

% === Format IEA data to standard format with ISO code ===
CO2_h_final_IEA(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','IEA_ISO'); 
error_index = strcmp(ISO_name_match(:, 2), 'ActiveX VT_ERROR: ');


ISO_name_match = ISO_name_match(~error_index, :);

CO2_hist_IEA=CO2_hist_IEA(2:end,:);
for i=1:length(CO2_hist_IEA)
    CO2_1=[CO2_hist_IEA(i,1),CO2_hist_IEA(i,12:22)];
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(1,2:end)];
    if length(aa)==13
    CO2_h_final_IEA=[CO2_h_final_IEA;aa];
        
    end
end
CO2_h_final_IEA(1,1:2)={'Country','ISO'};


% ===Combine own compiled emissions with IEA emissions===
CO2_h_final=[];
CO2_h_final_IEA=CO2_h_final_IEA(2:end,:);
for i=1:size(CO2_h_final_IEA,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1),:);
    else
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),CO2_h_final_IEA(i,2))==1),:);
    
    end
    
    CO2_h_final=[CO2_h_final;A];
end
[~,~,SSP_C]=xlsread('F:\research\code_and_data\data_generation\IAM_All_Country_Match(test) R5.xlsx','SSP-C'); 

[~,~,GDP_Country]=xlsread('F:\research\code_and_data\data_generation\SSP_scenarios.xlsx','SSP_scenarios (2)'); 
[~,~,POP_Country]=xlsread('F:\research\code_and_data\data_generation\SSP_scenarios.xlsx','SSP_scenarios (2)'); 


% === Extract and Filter GDP Data ===
coun_name_WDI=GDP_Country(find(ismember(GDP_Country(:,3),{'gdp'})==1),:);
filtered_GDP_Country = coun_name_WDI(find(ismember(coun_name_WDI(:,4),{'update'})==1),:);
filtered_GDP_Country = filtered_GDP_Country(find(ismember(filtered_GDP_Country(:,5),{'million 2005$PPP'})==1),:);
% === Extract and Filter POP Data ===
POP_Country=POP_Country(find(ismember(POP_Country(:,3),{'pop'})==1),:);
POP_Country_1=POP_Country(find(ismember(POP_Country(:,4),{'update'})==1),:);
filtered_POP_Country = POP_Country_1;
gdp_result = filtered_GDP_Country;

% === Interpolate Missing GDP Years===
GDP_result_cha=gdp_result;
for i=1:length(gdp_result)
xx=(2026:4:2030);
x = (2030:5:2100);
yy= cell2mat(gdp_result(i,52:53));
y = cell2mat(gdp_result(i,53:end));

xxi=2027:2029;
xi = 2031:1:2099;

yyi=pchip(xx,yy,xxi);
yi = pchip(x, y, xi);
yyc=[yy(1,1),yyi,yy(1,end)];
yc=[y(1,1),yi,y(1,end)];
GDP_result_cha(i,52:52+length(yyc)-1)=num2cell(yyc);
GDP_result_cha(i,52+length(yyc)-1:52+length(yyc)-1+length(yc)-1)=num2cell(yc);
end
% === Add Title Row and Format Final GDP Table ===
title=[{'iso3c'},{'SSP'},{'x'},{'version'},{'unit'},num2cell(1980:1:2100)];
GDP_result_cha=[title;GDP_result_cha];
GDP_result_cha=[GDP_result_cha(:,1:3),GDP_result_cha(:,5),GDP_result_cha(:,24:end)];


%=== Construct Matched CO2 and GDP Table ===
gdp_result_history=gdp_result(find(ismember(gdp_result(:,2),{'SSP1'})==1),:);
CO2_new_intensity=[];
for i=1:size(gdp_result_history,1)
    au=[gdp_result_history(i,1:2),gdp_result_history(i,36:46)];
    if ~isempty(find(ismember(CO2_h_final_own(:,2),gdp_result_history(i,1))==1))
    ui=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),gdp_result_history(i,1))==1),3:end);
    uii=[au(1,1),{'CO2'},ui];
    elseif ~isempty(find(ismember(CO2_h_final_IEA(:,2),gdp_result_history(i,1))==1))
      ui=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),gdp_result_history(i,1))==1),3:end);
    uii=[au(1,1),{'CO2'},ui];  
    else
        uii=[];
    end
    if isempty(uii)
        o=1;
    else
       uuii=[au;uii];
    CO2_new_intensity=[CO2_new_intensity;uuii];
    end
end
for i = 1:size(CO2_new_intensity, 1)
    for j = 1:size(CO2_new_intensity, 2)
        if isnan(CO2_new_intensity{i, j})
            CO2_new_intensity{i, j} = 0;
        end
    end
end

% === Compute Historical CO2 Intensity===
intensity_1=cell2mat(CO2_new_intensity(:,3:end));
intensity_final=[];
for i=1:size(intensity_1,1)./2
    intensity_2=intensity_1((i-1)*2+2,:)./intensity_1((i-1)*2+1,:);
    trans=[CO2_new_intensity((i-1)*2+2,1),num2cell(intensity_2)];
    intensity_final=[intensity_final;trans];
end
% ===Combine own compiled emissions with IEA emissions===
CO2_h_final=[];
for i=1:size(intensity_final,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),intensity_final(i,1))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),intensity_final(i,1))==1),:);
    elseif ~isempty(find(ismember(CO2_h_final_IEA(:,2),intensity_final(i,1))==1))
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),intensity_final(i,1))==1),:);
    else
        A=[];
    end
    
    CO2_h_final=[CO2_h_final;A];
end

yc=2100;
yh=2020;
r=log(0.01)/(yc-yh);
allmodel_matchhistory_CO2=[]; % Emission data that ultimately matches the total emissions of historical data

% ===IAMs results to individual countries===
for i=1:length(model_data)
    hu=model_data{i};
    
    [mapping_country_model, ~, country_match] = xlsread('F:\research\code_and_data\data_generation\IAM_All_Country_Match(test) R5.xlsx','AIM|CGE');
    mapping_country_model(isnan(mapping_country_model))=0;
    country_match=country_match(2:end,:);
    country_match(3:end-2,5:end)=num2cell(mapping_country_model);
    country_match(end-1:end,5:end)=num2cell(zeros(2, size(mapping_country_model, 2)));
    country_match=[country_match(:,1:4),country_match(:,find(ismember(country_match(1,:),hu(:,8))==1))];
    mapping_country_model=cell2mat(country_match(3:end,5:end));
    hangNo=find(sum(mapping_country_model,2)~=1&sum(mapping_country_model,2)~=0)+2;
    double_hang=country_match(hangNo,1:2);
    a=hu;
    a(strcmp(a, 'NA')) = {0};
    EIc_all_coun=[];% Intermediate results carbon intensity
    EIc_all_coun_final=[];% Intermediate results carbon intensity cali
    AIM_CGE_h_CO2=[];% Intermediate results aligned with historical emissions
    AIM_CGE_new_CO2=[];% Intermediate results not aligned with historical emissions
    AIM_CGE_new_GDP=[];%GDP
    for j=1:size(a,1)
        if sum(find(ismember(country_match(1,:),a(j,8))==1))~=0
            name_region=country_match(:,find(ismember(country_match(1,:),a(j,8))==1));
            if sum(cell2mat(name_region(3:end,1)))==1 %1 = country-level data, 0 = region data
             AIM_CGE_new_CO2=[AIM_CGE_new_CO2;a(j,:)];
             % ===country-level data matched to historical emissions===
             match_idx_1 = country_match(find(cell2mat(country_match(3:end,find(ismember(country_match(1,:),a(j,8))==1)))==1)+2,2);
             match_idx = find(strcmp(GDP_result_cha(:, 1), match_idx_1) & strcmp(GDP_result_cha(:, 2),SSP_C(strcmp(SSP_C(:,1),a(j,2)),3)));
             GDP_c=[a(j,1:8),'GDP','US',GDP_result_cha(match_idx,17:end)];
             AIM_CGE_new_GDP=[AIM_CGE_new_GDP; GDP_c];
             EIc=cell2mat(a(j,11:end))./cell2mat(GDP_c(:,11:end));
             combine=[a(j,1:10),num2cell(EIc)];
             EIc_all_coun=[EIc_all_coun;combine];

             % calculate annual growth rate relative to previous year
             b=a(j,:);
             growth_rate = diff(cell2mat(b(:,11:end)))./ cell2mat(b(:,11:end-1));
             b(2,11:end-1)=num2cell(growth_rate);
             % matched country-level data with historical data
             Country_CO2=CO2_h_final((find(ismember(CO2_h_final(:,2),match_idx_1)==1)),:);
             b(3,11:21)=Country_CO2(1,3:end);
             cc_1=cell2mat(b(2,11:end));
             cc_2=cell2mat(b(3,11:end));
             cc_1(2,1:11)=cc_2;
             for u=12:size(cc_1,2)+1
             cc_1(2,u) = cc_1(2,u-1) * (1+cc_1(1,u-1));
             end
             b(2:end,11:end)=num2cell(cc_1);
             b(3,1:10)=b(1,1:10);
             AIM_CGE_h_CO2=[AIM_CGE_h_CO2;b(3,:)];
             EIc_final=cell2mat(b(3,11:end))./cell2mat(GDP_c(:,11:end));
             combine=[a(j,1:10),num2cell(EIc_final)];
             EIc_all_coun_final=[EIc_all_coun_final;combine];
            else
            % matched regional data with historical data
             name_region=country_match(:,find(ismember(country_match(1,:),a(j,8))==1));
             name_V1=country_match(3:end,1:2);
             indivmodel_name_region=name_V1(find(cell2mat(name_region(3:end,1))==1),:);
             % if find(strcmp(name_region(1,1),country_match(1,:)))<9
             %  isMatched = ismember(indivmodel_name_region(:,:),double_hang(:,:));
             %  hui=reshape(indivmodel_name_region(isMatched), length(indivmodel_name_region(isMatched))./2,2);
             %  huii=country_match(find(ismember(country_match(:,2),hui(:,2))),:);
             %  [row, col] = find(cellfun(@(x) isnumeric(x) && x == 1, huii(:, 5:end)));
             %  col = col + 4;
             %  huiii=transpose(country_match(1,unique(col(col>9))));
             %  sumsum=sum(cell2mat(a(find(ismember(a(:,8),huiii(:,1))),11:end)),1);
             %  if isempty(sumsum)
             %    sumsum=zeros(1,size(a(j,11:end),2));
             %  end
             %  a(j,11:end)=num2cell(cell2mat(a(j,11:end))-sumsum);
             %  indivmodel_name_region(isMatched) = [];
             %  if any(isMatched(:) ~= 0)
             %  indivmodel_name_region=reshape(indivmodel_name_region, size(indivmodel_name_region,2)./2,2);
             %  end
             %  else
             %    indivmodel_name_region=indivmodel_name_region;
             %  end
              GDP_R=[];
              for h=1:size(indivmodel_name_region,1)
               match_idx = find(strcmp(GDP_result_cha(:,1), indivmodel_name_region(h,2)) & strcmp(GDP_result_cha(:, 2),SSP_C(strcmp(SSP_C(:,1),a(j,2)),3)));
               GDP_c=[a(j,1:7),GDP_result_cha(match_idx,1),'GDP','US',GDP_result_cha(match_idx,17:end)];
               GDP_R=[GDP_R;GDP_c];
              end 
            
              carbon_R=[];
              for h=1:size(indivmodel_name_region,1)
               match_idx = find(strcmp(CO2_h_final(:,2), indivmodel_name_region(h,2)));
               CO2_c=[a(j,1:7),CO2_h_final(match_idx,2),a(j,9:10),CO2_h_final(match_idx,3:end)];
               carbon_R=[carbon_R;CO2_c];
              end 
              carbon_R1=cell2mat(carbon_R(:,11:end));
              carbon_R1(isnan(carbon_R1))=0;
              carbon_R(:,11:end)=num2cell(carbon_R1);
             
              % calculate annual growth rate relative to previous year
              b=a(j,:);
              growth_rate = diff(cell2mat(b(:,11:end)))./ cell2mat(b(:,11:end-1));
              b(2,11:end-1)=num2cell(growth_rate);
              % Aggregate historical country-level data to region
              b(3,11:21)=num2cell(sum(cell2mat(carbon_R(:,11:end)),1));
              cc_1=cell2mat(b(2,11:end));
              cc_2=cell2mat(b(3,11:end));
              cc_1(2,1:11)=cc_2;
              for u=12:size(cc_1,2)+1
                cc_1(2,u) = cc_1(2,u-1) * (1+cc_1(1,u-1));
              end
              b(2:end,11:end)=num2cell(cc_1);
              b(3,1:10)=b(1,1:10);
              % record data
              GDP_region=sum(cell2mat(GDP_R(:,11:end)),1);
              GDP_region_final=[a(j,1:8),'GDP',a(j,10),num2cell(GDP_region)];
              intensity_r=num2cell(cell2mat(b(3,11:end))./cell2mat(GDP_region_final(:,11:end)));
              intensity_r_1=num2cell(cell2mat(a(j,11:end))./cell2mat(GDP_region_final(:,11:end)));
              AIM_CGE_new_CO2=[AIM_CGE_new_CO2;a(j,:)];
              AIM_CGE_h_CO2=[AIM_CGE_h_CO2;b(3,:)];
              AIM_CGE_new_GDP=[AIM_CGE_new_GDP;GDP_region_final];
              
              intensity_r_final=[a(j,1:10),intensity_r];
              intensity_r_1_final=[a(j,1:10),intensity_r_1];
              EIc_all_coun_final=[ EIc_all_coun_final;intensity_r_final];  
              EIc_all_coun=[ EIc_all_coun;intensity_r_1_final]; 
              % ===Convergence and target emission intensity===
              if isempty(find(cell2mat(b(3,11:end))<0))  
               for m=1:size(indivmodel_name_region,1)
                 % Calculate coefficient
                 count_intensity=intensity_final(find(ismember(intensity_final(:,1),indivmodel_name_region(m,2))==1),end);
                 ac=(cell2mat(intensity_r_final(:,end))-cell2mat(count_intensity))/(exp(r*yc)-exp(r*yh)); 
                 bc=(cell2mat(intensity_r_final(:,end))-0.01*cell2mat(count_intensity))/(1-0.01);
                 EIc_final=[];
                for y=2021:2100
                  % Calculate emission intensity in year y using formula
                  EIc=ac*exp(r*y)+bc;
                  EIc_final=[EIc_final,num2cell(EIc)];
                end  
                % Construction of the temporary emission intensity pathways
                match_idx = find(strcmp(intensity_final(:,1), indivmodel_name_region(m,2)));
                combine=[a(j,1:7),strjoin([name_region(1,1),indivmodel_name_region(m,2)], '_'),a(j,9:10), intensity_final(match_idx,2:12),EIc_final];
                EIc_all_coun=[ EIc_all_coun;combine];
                CO2_C=[combine(:,1:10),num2cell(cell2mat(GDP_R(m,11:end)).*cell2mat(combine(1,11:end)))];
                AIM_CGE_new_CO2=[AIM_CGE_new_CO2;CO2_C];
                AIM_CGE_new_GDP=[AIM_CGE_new_GDP;GDP_R(m,:)];
               end
              else
                yu=find(cell2mat(b(3,11:end))<0);
                yc_neg=2010+yu(1,1)-2;
                r_neg=log(0.01)/(yc_neg-yh);
                for m=1:size(indivmodel_name_region,1)
                    % Calculate coefficient
                    count_intensity=intensity_final(find(ismember(intensity_final(:,1),indivmodel_name_region(m,2))==1),end);
                    ac=(cell2mat(intensity_r_final(:,yu(1,1)+9))-cell2mat(count_intensity))/(exp(r_neg*yc_neg)-exp(r_neg*yh)); 
                    bc=(cell2mat(intensity_r_final(:,yu(1,1)+9))-0.01*cell2mat(count_intensity))/(1-0.01);
                    EIc_final=[];
                  for y=2021:yc_neg
                     % Calculate emission intensity in year y using formula
                     EIc=ac*exp(r_neg*y)+bc;
                     EIc_final=[EIc_final,num2cell(EIc)];
                  end
                % Construction of the temporary emission intensity pathways
                EIc_final(:,size(EIc_final,2)+1:80)=intensity_r_final(:,yu(1,1)+10:end);
                match_idx = find(strcmp(intensity_final(:,1), indivmodel_name_region(m,2)));
                combine=[a(j,1:7),strjoin([name_region(1,1),indivmodel_name_region(m,2)], '_'),a(j,9:10), intensity_final(match_idx,2:12),EIc_final];
                EIc_all_coun=[ EIc_all_coun;combine];
                CO2_C=[combine(:,1:10),num2cell(cell2mat(GDP_R(m,11:end)).*cell2mat(combine(1,11:end)))];
                AIM_CGE_new_CO2=[AIM_CGE_new_CO2;CO2_C];
                %AIM_CGE_new_GDP=[AIM_CGE_new_GDP;GDP_R(m,:)];
               end
              end
            % Emission pathways and scaling
            region_actual=AIM_CGE_h_CO2(end,:);
            region_CO2_temporary=AIM_CGE_new_CO2(end-size(GDP_R,1)+1:end,:);
            ratio=cell2mat(region_actual(:,11:end))./sum(cell2mat(region_CO2_temporary(:,11:end)),1);
            zuizhong=[region_CO2_temporary(:,1:10),num2cell(cell2mat(region_CO2_temporary(:,11:end)).*ratio)];
            AIM_CGE_h_CO2=[AIM_CGE_h_CO2;zuizhong];
            EI_zuizhong=[EIc_all_coun(end-size(GDP_R,1)+1:end,1:10),num2cell(cell2mat(zuizhong(:,11:end))./cell2mat(GDP_R(:,11:end)))];
            EIc_all_coun_final=[ EIc_all_coun_final;EI_zuizhong];
            end
        else
            bi=1;
        end
    end
        allmodel_matchhistory_CO2{i}=AIM_CGE_h_CO2; %final results
end
allmodel_matchhistory_CO2_final=vertcat(allmodel_matchhistory_CO2{:});%final results
xlswrite('F:\research\code_and_data\data_generation\final_result_scenarios.xlsx',allmodel_matchhistory_CO2_final,'CO2_match');

%% Inertia
clear all
cd 'F:\research\code_and_data\data_generation'
filename_1 = 'ISO_ar6_mainmodel.csv';  
CO2_ISO_total = readtable(filename_1, 'ReadVariableNames', true);
filename_2 = 'R5_ar6_mainmodel.csv';  
CO2_R5_total = readtable(filename_2, 'ReadVariableNames', true);
CO2_total_cha=[CO2_ISO_total;CO2_R5_total(2:end,:)];
% === Interpolation Module: Expand CO2 data from 5-year intervals to 1-year intervals ===
years = 2010:5:2100;
n_countries = height(CO2_total_cha);
result = array2table([years'], 'VariableNames', {'Year'});
CO2_total_cha_final=CO2_total_cha(:,1:10);
new_years = 2010:2100;
new_year_names = strcat("Y", string(new_years));
CO2_total_cha_final(:,11:101) = array2table(nan(height(CO2_total_cha), numel(new_years)));  
CO2_total_cha_final.Properties.VariableNames(11:101) = new_year_names;

for i = 1:n_countries
   
    
 
    row_data = table2array(CO2_total_cha(i, 11:end)); 
    year_idx = years(~isnan(row_data)); 
    values = row_data(~isnan(row_data));
    
   
    
    
    % === interp ===
    interp_years = 2010:1:max(year_idx);
    interp_values = pchip(year_idx, values, interp_years);
    CO2_total_cha_final(i,11:end)=array2table(interp_values);
   
    
end
% === Convert Table to Cell Format ===
% Extract variable names as a cell row
varNames = CO2_total_cha_final.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "Y");
dataCell = table2cell(CO2_total_cha_final);
CO2_total_cha_final_cell = [varNames; dataCell];
CO2_total_cha=CO2_total_cha_final_cell;


varNames = CO2_ISO_total.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "x");
dataCell = table2cell(CO2_ISO_total);
CO2_ISO_total = [varNames; dataCell];


varNames = CO2_R5_total.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "x");
dataCell = table2cell(CO2_R5_total);
CO2_R5_total = [varNames; dataCell];
clearvars -except CO2_total_cha CO2_ISO_total CO2_R5_total


models_name = CO2_R5_total(2:end,1);
models_name_1 = CO2_total_cha(:,1);

unique_models_name = unique(models_name);

% === Group Data by Model ===
% For each unique model, extract all matching rows from CO2_total_cha and store in cell array
model_data = cell(size(unique_models_name));
model_data_name=[];
for i = 1:numel(unique_models_name)
    model = unique_models_name(i);
    model_data{i}=CO2_total_cha(find(ismember(models_name_1,model)==1),:);
    model_data_name=[model_data_name;model];
end

[~,~,CO2_hist_own]=xlsread('F:\research\code_and_data\data_generation\Comparison of emissions for 60 developing countries with IEA.xlsx'); 
[~,~,CO2_hist_IEA]=xlsread('F:\research\code_and_data\data_generation\institution_comparison.xlsx','IEA_MT'); 
CO2_hist_own_summedMatrix=[CO2_hist_own(:,1:3);CO2_hist_own(2:end,5:7)];
CO2_hist_IEA=transpose(CO2_hist_IEA);

% === Format historical own data to standard format with ISO code ===
CO2_h_final_own(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','own_ISO'); 
CO2_coun_unique = unique(CO2_hist_own_summedMatrix(:,1));
for i=1:length(CO2_coun_unique)
    CO2_1=CO2_hist_own_summedMatrix(find(ismember(CO2_hist_own_summedMatrix(:,1),CO2_coun_unique(i,:))==1),:);
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(:,3)'];
    if length(aa)==13
    CO2_h_final_own=[CO2_h_final_own;aa];
        
    end
end
CO2_h_final_own(1,1:2)={'Country','ISO'};
CO2_h_final_own=CO2_h_final_own(:,1:13);

% === Format IEA data to standard format with ISO code ===
CO2_h_final_IEA(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','IEA_ISO'); 
error_index = strcmp(ISO_name_match(:, 2), 'ActiveX VT_ERROR: ');


ISO_name_match = ISO_name_match(~error_index, :);

CO2_hist_IEA=CO2_hist_IEA(2:end,:);
for i=1:length(CO2_hist_IEA)
    CO2_1=[CO2_hist_IEA(i,1),CO2_hist_IEA(i,12:22)];
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(1,2:end)];
    if length(aa)==13
    CO2_h_final_IEA=[CO2_h_final_IEA;aa];
        
    end
end
CO2_h_final_IEA(1,1:2)={'Country','ISO'};


% ===Combine own compiled emissions with IEA emissions===
CO2_h_final=[];
CO2_h_final_IEA=CO2_h_final_IEA(2:end,:);
for i=1:size(CO2_h_final_IEA,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1),:);
    else
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),CO2_h_final_IEA(i,2))==1),:);
    
    end
    
    CO2_h_final=[CO2_h_final;A];
end
[~,~,SSP_C]=xlsread('F:\research\code_and_data\data_generation\IAM_All_Country_Match(test) R5.xlsx','SSP-C'); 

[~,~,GDP_Country]=xlsread('F:\research\code_and_data\data_generation\SSP_scenarios.xlsx','SSP_scenarios (2)'); 
[~,~,POP_Country]=xlsread('F:\research\code_and_data\data_generation\SSP_scenarios.xlsx','SSP_scenarios (2)'); 


% === Extract and Filter GDP Data ===
coun_name_WDI=GDP_Country(find(ismember(GDP_Country(:,3),{'gdp'})==1),:);
filtered_GDP_Country = coun_name_WDI(find(ismember(coun_name_WDI(:,4),{'update'})==1),:);
filtered_GDP_Country = filtered_GDP_Country(find(ismember(filtered_GDP_Country(:,5),{'million 2005$PPP'})==1),:);
% === Extract and Filter POP Data ===
POP_Country=POP_Country(find(ismember(POP_Country(:,3),{'pop'})==1),:);
POP_Country_1=POP_Country(find(ismember(POP_Country(:,4),{'update'})==1),:);
filtered_POP_Country = POP_Country_1;
gdp_result = filtered_GDP_Country;

% === Interpolate Missing GDP Years===
GDP_result_cha=gdp_result;
for i=1:length(gdp_result)
xx=(2026:4:2030);
x = (2030:5:2100);
yy= cell2mat(gdp_result(i,52:53));
y = cell2mat(gdp_result(i,53:end));

xxi=2027:2029;
xi = 2031:1:2099;

yyi=pchip(xx,yy,xxi);
yi = pchip(x, y, xi);
yyc=[yy(1,1),yyi,yy(1,end)];
yc=[y(1,1),yi,y(1,end)];
GDP_result_cha(i,52:52+length(yyc)-1)=num2cell(yyc);
GDP_result_cha(i,52+length(yyc)-1:52+length(yyc)-1+length(yc)-1)=num2cell(yc);
end
% === Add Title Row and Format Final GDP Table ===
title=[{'iso3c'},{'SSP'},{'x'},{'version'},{'unit'},num2cell(1980:1:2100)];
GDP_result_cha=[title;GDP_result_cha];
GDP_result_cha=[GDP_result_cha(:,1:3),GDP_result_cha(:,5),GDP_result_cha(:,24:end)];


%=== Construct Matched CO2 and GDP Table ===
gdp_result_history=gdp_result(find(ismember(gdp_result(:,2),{'SSP1'})==1),:);
CO2_new_intensity=[];
for i=1:size(gdp_result_history,1)
    au=[gdp_result_history(i,1:2),gdp_result_history(i,36:46)];
    if ~isempty(find(ismember(CO2_h_final_own(:,2),gdp_result_history(i,1))==1))
    ui=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),gdp_result_history(i,1))==1),3:end);
    uii=[au(1,1),{'CO2'},ui];
    elseif ~isempty(find(ismember(CO2_h_final_IEA(:,2),gdp_result_history(i,1))==1))
      ui=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),gdp_result_history(i,1))==1),3:end);
    uii=[au(1,1),{'CO2'},ui];  
    else
        uii=[];
    end
    if isempty(uii)
        o=1;
    else
       uuii=[au;uii];
    CO2_new_intensity=[CO2_new_intensity;uuii];
    end
end
for i = 1:size(CO2_new_intensity, 1)
    for j = 1:size(CO2_new_intensity, 2)
        if isnan(CO2_new_intensity{i, j})
            CO2_new_intensity{i, j} = 0;
        end
    end
end

% === Compute Historical CO2 Intensity===
intensity_1=cell2mat(CO2_new_intensity(:,3:end));
intensity_final=[];
for i=1:size(intensity_1,1)./2
    intensity_2=intensity_1((i-1)*2+2,:)./intensity_1((i-1)*2+1,:);
    trans=[CO2_new_intensity((i-1)*2+2,1),num2cell(intensity_2)];
    intensity_final=[intensity_final;trans];
end
% ===Combine own compiled emissions with IEA emissions===
CO2_h_final=[];
for i=1:size(intensity_final,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),intensity_final(i,1))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),intensity_final(i,1))==1),:);
    elseif ~isempty(find(ismember(CO2_h_final_IEA(:,2),intensity_final(i,1))==1))
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),intensity_final(i,1))==1),:);
    else
        A=[];
    end
    
    CO2_h_final=[CO2_h_final;A];
end
clearvars -except CO2_h_final intensity_final GDP_result_cha

GDP_result_cha_SSP2=GDP_result_cha(find(ismember(GDP_result_cha(:,2),{'SSP2'})==1),:);
GDP_history_high=(cell2mat(GDP_result_cha_SSP2(:,26))./cell2mat(GDP_result_cha_SSP2(:,22))).^(1/4)-1;
GDP_history_high=repmat(GDP_history_high,1,30);
GDP_2020=cell2mat(GDP_result_cha_SSP2(:,27));
GDP_2020=repmat(GDP_2020,1,30);
years = (2021:2050); 
num_years = length(years); 
GDP_HIS_H = [GDP_2020.* cumprod(1 +GDP_history_high,2)];
GDP_result_cha_SSP2(:,28:57)=num2cell(GDP_HIS_H);
% === Estimate Emissions using Average Intensity ===
BAU_average_final=[];
for i=1:size(CO2_h_final,1)
    a=GDP_result_cha_SSP2(find(ismember(GDP_result_cha_SSP2(:,1),CO2_h_final(i,2))==1),:);
    average=sum(cell2mat(intensity_final(i,7:end-1)),2)./5;
    BAU_average=[CO2_h_final(i,:),num2cell(average.*cell2mat(a(:,28:57)))];
    BAU_average_final=[BAU_average_final;BAU_average];
     
end


% === Estimate Emissions using Annual Intensity Decline ===
BAU_annual_average_final=[];
for i=1:size(CO2_h_final,1)
    a=GDP_result_cha_SSP2(find(ismember(GDP_result_cha_SSP2(:,1),CO2_h_final(i,2))==1),:);
 
    % Calculate annual intensity growth rate
    annual_growth_rate = (cell2mat(intensity_final(i,end-1)) /cell2mat(intensity_final(i,7)))^(1 / 4) - 1;
    if annual_growth_rate<=0
    annual_growth_factor = (1 + annual_growth_rate);
    years = 2050 - 2020+1; 
    intensity_t = zeros(1, years);
    intensity_t(1) = cell2mat(intensity_final(i,end)); % 2020

    for j = 2:years
    intensity_t(j) = intensity_t(j-1) * annual_growth_factor;
    end
   
    BAU_annual_average=[CO2_h_final(i,:),num2cell(intensity_t(:,2:end).*cell2mat(a(:,28:57)))];
    BAU_annual_average_final=[BAU_annual_average_final;BAU_annual_average];  
    else
        % If intensity is increasing, fallback to average method
        BAU_annual_average_final=[BAU_annual_average_final;BAU_average_final(i,:)];
    end
end


xlswrite('F:\research\code_and_data\data_generation\BAU.xlsx',BAU_annual_average_final);

%% excellence scenarios
clear all
cd 'F:\research\code_and_data\data_generation'
filename_1 = 'ISO_ar6_mainmodel.csv';  
CO2_ISO_total = readtable(filename_1, 'ReadVariableNames', true);
filename_2 = 'R5_ar6_mainmodel.csv';  
CO2_R5_total = readtable(filename_2, 'ReadVariableNames', true);
CO2_total_cha=[CO2_ISO_total;CO2_R5_total(2:end,:)];
% === Interpolation Module: Expand CO2 data from 5-year intervals to 1-year intervals ===
years = 2010:5:2100;
n_countries = height(CO2_total_cha);
result = array2table([years'], 'VariableNames', {'Year'});
CO2_total_cha_final=CO2_total_cha(:,1:10);
new_years = 2010:2100;
new_year_names = strcat("Y", string(new_years));
CO2_total_cha_final(:,11:101) = array2table(nan(height(CO2_total_cha), numel(new_years)));  
CO2_total_cha_final.Properties.VariableNames(11:101) = new_year_names;

for i = 1:n_countries
   
    
 
    row_data = table2array(CO2_total_cha(i, 11:end)); 
    year_idx = years(~isnan(row_data)); 
    values = row_data(~isnan(row_data));
    
   
    
    
    % === interp ===
    interp_years = 2010:1:max(year_idx);
    interp_values = pchip(year_idx, values, interp_years);
    CO2_total_cha_final(i,11:end)=array2table(interp_values);
   
    
end
% === Convert Table to Cell Format ===
% Extract variable names as a cell row
varNames = CO2_total_cha_final.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "Y");
dataCell = table2cell(CO2_total_cha_final);
CO2_total_cha_final_cell = [varNames; dataCell];
CO2_total_cha=CO2_total_cha_final_cell;


varNames = CO2_ISO_total.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "x");
dataCell = table2cell(CO2_ISO_total);
CO2_ISO_total = [varNames; dataCell];


varNames = CO2_R5_total.Properties.VariableNames;
varNames(11:end) = extractAfter(varNames(11:end), "x");
dataCell = table2cell(CO2_R5_total);
CO2_R5_total = [varNames; dataCell];
clearvars -except CO2_total_cha CO2_ISO_total CO2_R5_total


models_name = CO2_R5_total(2:end,1);
models_name_1 = CO2_total_cha(:,1);

unique_models_name = unique(models_name);

% === Group Data by Model ===
% For each unique model, extract all matching rows from CO2_total_cha and store in cell array
model_data = cell(size(unique_models_name));
model_data_name=[];
for i = 1:numel(unique_models_name)
    model = unique_models_name(i);
    model_data{i}=CO2_total_cha(find(ismember(models_name_1,model)==1),:);
    model_data_name=[model_data_name;model];
end

[~,~,CO2_hist_own]=xlsread('F:\research\code_and_data\data_generation\Comparison of emissions for 60 developing countries with IEA.xlsx'); 
[~,~,CO2_hist_IEA]=xlsread('F:\research\code_and_data\data_generation\institution_comparison.xlsx','IEA_MT'); 
CO2_hist_own_summedMatrix=[CO2_hist_own(:,1:3);CO2_hist_own(2:end,5:7)];
CO2_hist_IEA=transpose(CO2_hist_IEA);

% === Format historical own data to standard format with ISO code ===
CO2_h_final_own(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','own_ISO'); 
CO2_coun_unique = unique(CO2_hist_own_summedMatrix(:,1));
for i=1:length(CO2_coun_unique)
    CO2_1=CO2_hist_own_summedMatrix(find(ismember(CO2_hist_own_summedMatrix(:,1),CO2_coun_unique(i,:))==1),:);
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(:,3)'];
    if length(aa)==13
    CO2_h_final_own=[CO2_h_final_own;aa];
        
    end
end
CO2_h_final_own(1,1:2)={'Country','ISO'};
CO2_h_final_own=CO2_h_final_own(:,1:13);

% === Format IEA data to standard format with ISO code ===
CO2_h_final_IEA(1,3:13)=num2cell(2010:1:2020);
[~,~,ISO_name_match]=xlsread('F:\research\code_and_data\data_generation\ISO3_Codes_Country_Names.xlsx','IEA_ISO'); 
error_index = strcmp(ISO_name_match(:, 2), 'ActiveX VT_ERROR: ');


ISO_name_match = ISO_name_match(~error_index, :);

CO2_hist_IEA=CO2_hist_IEA(2:end,:);
for i=1:length(CO2_hist_IEA)
    CO2_1=[CO2_hist_IEA(i,1),CO2_hist_IEA(i,12:22)];
    aa=[CO2_1(1,1),ISO_name_match(find(ismember(ISO_name_match(:,1),CO2_1(1,1))==1),2),CO2_1(1,2:end)];
    if length(aa)==13
    CO2_h_final_IEA=[CO2_h_final_IEA;aa];
        
    end
end
CO2_h_final_IEA(1,1:2)={'Country','ISO'};


% ===Combine own compiled emissions with IEA emissions===
CO2_h_final=[];
CO2_h_final_IEA=CO2_h_final_IEA(2:end,:);
for i=1:size(CO2_h_final_IEA,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),CO2_h_final_IEA(i,2))==1),:);
    else
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),CO2_h_final_IEA(i,2))==1),:);
    
    end
    
    CO2_h_final=[CO2_h_final;A];
end
[~,~,SSP_C]=xlsread('F:\research\code_and_data\data_generation\IAM_All_Country_Match(test) R5.xlsx','SSP-C'); 

[~,~,GDP_Country]=xlsread('F:\research\code_and_data\data_generation\SSP_scenarios.xlsx','SSP_scenarios (2)'); 
[~,~,POP_Country]=xlsread('F:\research\code_and_data\data_generation\SSP_scenarios.xlsx','SSP_scenarios (2)'); 
[~,~,gdppc_Country]=xlsread('F:\research\code_and_data\data_generation\SSP_scenarios.xlsx','SSP_scenarios (2)'); 


% === Extract and Filter GDP Data ===
coun_name_WDI=GDP_Country(find(ismember(GDP_Country(:,3),{'gdp'})==1),:);
filtered_GDP_Country = coun_name_WDI(find(ismember(coun_name_WDI(:,4),{'update'})==1),:);
filtered_GDP_Country = filtered_GDP_Country(find(ismember(filtered_GDP_Country(:,5),{'million 2005$PPP'})==1),:);
% === Extract and Filter POP Data ===
POP_Country=POP_Country(find(ismember(POP_Country(:,3),{'pop'})==1),:);
POP_Country_1=POP_Country(find(ismember(POP_Country(:,4),{'update'})==1),:);
filtered_POP_Country = POP_Country_1;
filtered_POP_Country=[filtered_POP_Country(:,1:3),filtered_POP_Country(:,5:end)];

gdp_result = filtered_GDP_Country;

% === Extract and Filter GDP Per Capita Data ===
gdppc_Country=gdppc_Country(find(ismember(gdppc_Country(:,3),{'gdppc'})==1),:);
gdppc_Country_1=gdppc_Country(find(ismember(gdppc_Country(:,4),{'update'})==1),:);
gdppc_Country_1 = gdppc_Country_1(find(ismember(gdppc_Country_1(:,5),{'2005$PPP'})==1),:);
gdppc_result = gdppc_Country_1;

% === Interpolate Missing GDP Per Capita Years===
GDPPC_result_cha=gdppc_result;
for i=1:length(gdppc_result)
xx=(2026:4:2030);
x = (2030:5:2100);
yy= cell2mat(gdppc_result(i,52:53));
y = cell2mat(gdppc_result(i,53:end));

xxi=2027:2029;
xi = 2031:1:2099;

yyi=pchip(xx,yy,xxi);
yi = pchip(x, y, xi);
yyc=[yy(1,1),yyi,yy(1,end)];
yc=[y(1,1),yi,y(1,end)];
GDPPC_result_cha(i,52:52+length(yyc)-1)=num2cell(yyc);
GDPPC_result_cha(i,52+length(yyc)-1:52+length(yyc)-1+length(yc)-1)=num2cell(yc);
end
% === Add Title Row and Format Final GDP Per Capita Table ===
title=[{'iso3c'},{'SSP'},{'x'},{'version'},{'unit'},num2cell(1980:1:2100)];
GDPPC_result_cha=[title;GDPPC_result_cha];
GDPPC_result_cha=[GDPPC_result_cha(:,1:3),GDPPC_result_cha(:,5),GDPPC_result_cha(:,24:end)];


% === Interpolate Missing GDP Years===
GDP_result_cha=gdp_result;
for i=1:length(gdp_result)
xx=(2026:4:2030);
x = (2030:5:2100);
yy= cell2mat(gdp_result(i,52:53));
y = cell2mat(gdp_result(i,53:end));

xxi=2027:2029;
xi = 2031:1:2099;

yyi=pchip(xx,yy,xxi);
yi = pchip(x, y, xi);
yyc=[yy(1,1),yyi,yy(1,end)];
yc=[y(1,1),yi,y(1,end)];
GDP_result_cha(i,52:52+length(yyc)-1)=num2cell(yyc);
GDP_result_cha(i,52+length(yyc)-1:52+length(yyc)-1+length(yc)-1)=num2cell(yc);
end
% === Add Title Row and Format Final GDP Table ===
title=[{'iso3c'},{'SSP'},{'x'},{'version'},{'unit'},num2cell(1980:1:2100)];
GDP_result_cha=[title;GDP_result_cha];
GDP_result_cha=[GDP_result_cha(:,1:3),GDP_result_cha(:,5),GDP_result_cha(:,24:end)];

%=== Construct Matched CO2 and GDP Table ===
gdp_result_history=gdp_result(find(ismember(gdp_result(:,2),{'SSP1'})==1),:);
CO2_new_intensity=[];
for i=1:size(gdp_result_history,1)
    au=[gdp_result_history(i,1:2),gdp_result_history(i,36:46)];
    if ~isempty(find(ismember(CO2_h_final_own(:,2),gdp_result_history(i,1))==1))
    ui=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),gdp_result_history(i,1))==1),3:end);
    uii=[au(1,1),{'CO2'},ui];
    elseif ~isempty(find(ismember(CO2_h_final_IEA(:,2),gdp_result_history(i,1))==1))
      ui=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),gdp_result_history(i,1))==1),3:end);
    uii=[au(1,1),{'CO2'},ui];  
    else
        uii=[];
    end
    if isempty(uii)
        o=1;
    else
       uuii=[au;uii];
    CO2_new_intensity=[CO2_new_intensity;uuii];
    end
end
for i = 1:size(CO2_new_intensity, 1)
    for j = 1:size(CO2_new_intensity, 2)
        if isnan(CO2_new_intensity{i, j})
            CO2_new_intensity{i, j} = 0;
        end
    end
end

% === Compute Historical CO2 Intensity===
intensity_1=cell2mat(CO2_new_intensity(:,3:end));
intensity_final=[];
for i=1:size(intensity_1,1)./2
    intensity_2=intensity_1((i-1)*2+2,:)./intensity_1((i-1)*2+1,:);
    trans=[CO2_new_intensity((i-1)*2+2,1),num2cell(intensity_2)];
    intensity_final=[intensity_final;trans];
end
% ===Combine own compiled emissions with IEA emissions===
CO2_h_final=[];
for i=1:size(intensity_final,1)
    if ~isempty(find(ismember(CO2_h_final_own(:,2),intensity_final(i,1))==1))
    A=CO2_h_final_own(find(ismember(CO2_h_final_own(:,2),intensity_final(i,1))==1),:);
    elseif ~isempty(find(ismember(CO2_h_final_IEA(:,2),intensity_final(i,1))==1))
     A=CO2_h_final_IEA(find(ismember(CO2_h_final_IEA(:,2),intensity_final(i,1))==1),:);
    else
        A=[];
    end
    
    CO2_h_final=[CO2_h_final;A];
end
clearvars -except CO2_h_final intensity_final GDP_result_cha GDPPC_result_cha
GDPPC_result_cha_SSP2=GDPPC_result_cha(find(ismember(GDPPC_result_cha(:,2),{'SSP2'})==1),:);
gdppc_final = cell(size(intensity_final));


for i = 1:size(intensity_final, 1)
    country = intensity_final{i, 1};
    
    
    match_idx = find(strcmp(GDPPC_result_cha_SSP2(:, 1), country));
    gdppc = GDPPC_result_cha_SSP2(match_idx, 17:27);
    
    
    gdppc_final{i, 1} = country;
    gdppc_final(i, 2:end) = gdppc;
end
clearvars -except CO2_h_final intensity_final GDP_result_cha GDPPC_result_cha gdppc_final

% ===Data Preparation===
intensity_change_1=(cell2mat(intensity_final(:,7))./cell2mat(intensity_final(:,2))).^(1/5)-1;
intensity_change_2=(cell2mat(intensity_final(:,11))./cell2mat(intensity_final(:,7))).^(1/4)-1;
intensity_change_3=[intensity_final(:,1),num2cell(intensity_change_1),num2cell(intensity_change_2)];
gdppc_change_1=(cell2mat(gdppc_final(:,7))./cell2mat(gdppc_final(:,2))).^(1/5)-1;
gdppc_change_2=(cell2mat(gdppc_final(:,11))./cell2mat(gdppc_final(:,7))).^(1/4)-1;
gdppc_change_3=[gdppc_final(:,1),num2cell(gdppc_change_1),num2cell(gdppc_change_2)];
index=cell2mat(gdppc_change_3(:,2:3))./cell2mat(intensity_change_3(:,2:3));
total=[gdppc_change_3,intensity_change_3(:,2:3),num2cell(index)];
clearvars -except CO2_h_final intensity_final GDP_result_cha GDPPC_result_cha gdppc_final total

GDP_result_cha_SSP2=GDP_result_cha(find(ismember(GDP_result_cha(:,2),{'SSP2'})==1),:);
GDP_result_cha_SSP2_v2 = cell(size(intensity_final));


for i = 1:size(intensity_final, 1)
    country = intensity_final{i, 1};
    
   
    match_idx = find(strcmp(GDP_result_cha_SSP2(:, 1), country));
    gdppc = GDP_result_cha_SSP2(match_idx, 17:27);
    
    
    GDP_result_cha_SSP2_v2 {i, 1} = country;
    GDP_result_cha_SSP2_v2 (i, 2:end) = gdppc;
end
clearvars -except CO2_h_final intensity_final GDP_result_cha GDPPC_result_cha gdppc_final total GDP_result_cha_SSP2_v2
GDP_history_high=(cell2mat(GDP_result_cha_SSP2_v2(:,11))./cell2mat(GDP_result_cha_SSP2_v2(:,7))).^(1/4)-1;%GDP annual average change rate
GDP_change_rate=[gdppc_final(:,1),num2cell(GDP_history_high)];
GDP_history_high=repmat(GDP_history_high,1,30);
GDP_2020=cell2mat(GDP_result_cha_SSP2_v2(:,end));
GDP_2020=repmat(GDP_2020,1,30);
years = (2021:2050); 
num_years = length(years); 
GDP_HIS_H = [GDP_2020.* cumprod(1 +GDP_history_high,2)];
GDP_result_cha_SSP2_v2(:,13:42)=num2cell(GDP_HIS_H);
clearvars -except CO2_h_final intensity_final GDP_result_cha GDPPC_result_cha gdppc_final total GDP_result_cha_SSP2_v2 GDP_change_rate
GDPPC_history_high=repmat(total(:,3),1,30);
GDPPC_2020=cell2mat(gdppc_final(:,end));
GDPPC_2020=repmat(GDPPC_2020,1,30);
years = (2021:2050); 
num_years = length(years); 
GDPPC_HIS_H = [GDPPC_2020.* cumprod(1 +cell2mat(GDPPC_history_high),2)];
gdppc_final(:,13:42)=num2cell(GDPPC_HIS_H);
clearvars -except CO2_h_final intensity_final GDP_result_cha GDPPC_result_cha gdppc_final total GDP_result_cha_SSP2_v2 GDP_change_rate
GDP_result_cha_SSP2_v2=GDP_result_cha_SSP2_v2(:,1:12);
gdppc_final=gdppc_final(:,1:12);
a=[gdppc_final(:,1),gdppc_final(:,12)];

% ===Country Ordering===
sorted_data = sortrows(a, 2);
total_rows = size(sorted_data, 1); 

intensity_final_final=[];
gdppc_final_final=[];
GDP_final_final=[];
CO2_final_final=[];
for i=1:size(gdppc_final,1)
    country=gdppc_final(i,1);
    target_2020=gdppc_final(i,end);
    values = cell2mat(sorted_data(:,2)); 
    % ===Peer Selection===
    [~, closest_idx] = min(abs(values - cell2mat(target_2020))); 
    if closest_idx <= 5
        start_idx=1;
        end_idx=11;
    elseif closest_idx > total_rows - 5   
        start_idx=total_rows-10;
        end_idx=total_rows;
    else
    start_idx = closest_idx - 5; 
    end_idx = closest_idx + 5; 
    end
     countries=[sorted_data(start_idx:closest_idx-1,:);sorted_data(closest_idx+1:end_idx,:)];
     %change rates of GDP, intensity, GDP per capita
    for j=1:size(countries,1)
        countries(j,3)=total(find(ismember(total(:,1),countries(j,1))==1),3);
        countries(j,4)=total(find(ismember(total(:,1),countries(j,1))==1),5);
        countries(j,5)=GDP_change_rate(find(ismember(GDP_change_rate(:,1),countries(j,1))==1),2);
    end
    % ===Peer filtering by performance===
     country_own=[total(i,1),total(i,3),total(i,5),GDP_change_rate(i,2)];
     countries=sortrows(countries, 4);
     countries = countries(cell2mat(countries(:, 5)) > 0, :);
     countries = countries(cell2mat(countries(:, 4)) < cell2mat(country_own(:,3)), :);
     countries = countries(cell2mat(countries(:, 4)) < 0, :);
     countries(:,6)=num2cell(cell2mat(countries(:,5))./abs(cell2mat(countries(:,4))));
     countries=sortrows(countries, 6, 'descend');
     if isempty(countries)
         countries=[country_own(:,1),{0},country_own(:,2:end)];
     end
     average_intensity=repmat(cell2mat(countries(:,4)),1,5);
     average_GDP=repmat(cell2mat(countries(:,5)),1,5);
     average_gdppc=repmat(cell2mat(countries(:,2)),1,5);

     intensity= num2cell(cell2mat(intensity_final(i,12)).* cumprod(1 +average_intensity,2));
     combine=[repmat(intensity_final(i,:) ,size(intensity,1),1),intensity];
     intensity_final_final=[intensity_final_final;combine];
     gdppc=num2cell(cell2mat(gdppc_final(i,12)).* cumprod(1 +average_gdppc,2));
     combine_gdppc=[repmat(gdppc_final(i,:) ,size(gdppc,1),1),gdppc];
     gdppc_final_final=[gdppc_final_final;combine_gdppc];
    
     GDP=num2cell(cell2mat(GDP_result_cha_SSP2_v2(i,12)).* cumprod(1 +average_GDP,2));
     combine_GDP=[repmat(GDP_result_cha_SSP2_v2(i,:) ,size(GDP,1),1),GDP];
     GDP_final_final=[GDP_final_final;combine_GDP];
     CO2=num2cell(cell2mat(GDP_final_final((size(GDP_final_final,1)-size(countries,1)+1):end,13:17)).*cell2mat(intensity_final_final((size(GDP_final_final,1)-size(countries,1)+1):end,13:17)));
     combine_CO2=[repmat(CO2_h_final(i,:) ,size(CO2,1),1),CO2];
     CO2_final_final=[CO2_final_final;combine_CO2];
  
end
% ===Iterate Every Five Years===
intensity_final_2030=[];
gdppc_final_2030=[];
GDP_final_2030=[];
CO2_final_2030=[];
for i=1:size(gdppc_final_final,1)
    country=gdppc_final_final(i,1);
    target_2020=gdppc_final_final(i,end);
    values = cell2mat(sorted_data(:,2)); 
    % ===Peer Selection===
    [~, closest_idx] = min(abs(values - cell2mat(target_2020))); 
    if closest_idx <= 5
        start_idx=1;
        end_idx=11;
    elseif closest_idx > total_rows - 5   
        start_idx=total_rows-10;
        end_idx=total_rows;
    else
    start_idx = closest_idx - 5; 
    end_idx = closest_idx + 5; 
    end
     countries=[sorted_data(start_idx:closest_idx-1,:);sorted_data(closest_idx+1:end_idx,:)];
    
    for j=1:size(countries,1)
        countries(j,3)=total(find(ismember(total(:,1),countries(j,1))==1),3);
        countries(j,4)=total(find(ismember(total(:,1),countries(j,1))==1),5);
        countries(j,5)=GDP_change_rate(find(ismember(GDP_change_rate(:,1),countries(j,1))==1),2);
    end
    % ===Peer filtering by performance===
         country_own=[total(find(ismember(total(:,1),country)==1),1),total(find(ismember(total(:,1),country)==1),3),total(find(ismember(total(:,1),country)==1),5),GDP_change_rate(find(ismember(GDP_change_rate(:,1),country)==1),2)];
         countries=sortrows(countries, 4);
         countries = countries(cell2mat(countries(:, 5)) > 0, :);
         countries = countries(cell2mat(countries(:, 4)) < cell2mat(country_own(:,3)), :);
         countries = countries(cell2mat(countries(:, 4)) < 0, :);
         countries(:,6)=num2cell(cell2mat(countries(:,5))./abs(cell2mat(countries(:,4))));
         countries=sortrows(countries, 6, 'descend');
     if isempty(countries)
         countries=[country_own(:,1),{0},country_own(:,2:end)];
     end

      average_intensity=repmat(cell2mat(countries(:,4)),1,5);
      average_GDP=repmat(cell2mat(countries(:,5)),1,5);
      average_gdppc=repmat(cell2mat(countries(:,2)),1,5);

      intensity= num2cell(cell2mat(intensity_final_final(i,end)).* cumprod(1 +average_intensity,2));
      combine=[repmat(intensity_final_final(i,:) ,size(intensity,1),1),intensity];
      intensity_final_2030=[intensity_final_2030;combine];
      gdppc=num2cell(cell2mat(gdppc_final_final(i,17)).* cumprod(1 +average_gdppc,2));
      combine_gdppc=[repmat(gdppc_final_final(i,:) ,size(gdppc,1),1),gdppc];
      gdppc_final_2030=[gdppc_final_2030;combine_gdppc];
    
      GDP=num2cell(cell2mat(GDP_final_final(i,17)).* cumprod(1 +average_GDP,2));
      combine_GDP=[repmat(GDP_final_final(i,:) ,size(GDP,1),1),GDP];
      GDP_final_2030=[GDP_final_2030;combine_GDP];
      CO2=num2cell(cell2mat(GDP_final_2030((size(GDP_final_2030,1)-size(countries,1)+1):end,18:22)).*cell2mat(intensity_final_2030((size(GDP_final_2030,1)-size(countries,1)+1):end,18:22)));
      combine_CO2=[repmat(CO2_final_final(i,:) ,size(CO2,1),1),CO2];
      CO2_final_2030=[CO2_final_2030;combine_CO2];

end
intensity_final_2035=[];
gdppc_final_2035=[];
GDP_final_2035=[];
CO2_final_2035=[];
for i=1:size(gdppc_final_2030,1)
    country=gdppc_final_2030(i,1);
    target_2020=gdppc_final_2030(i,end);
    values = cell2mat(sorted_data(:,2)); 
    % ===Peer Selection===
    [~, closest_idx] = min(abs(values - cell2mat(target_2020))); 
    if closest_idx <= 5
        start_idx=1;
        end_idx=11;
    elseif closest_idx > total_rows - 5   
        start_idx=total_rows-10;
        end_idx=total_rows;
    else
    start_idx = closest_idx - 5; 
    end_idx = closest_idx + 5; 
    end
     countries=[sorted_data(start_idx:closest_idx-1,:);sorted_data(closest_idx+1:end_idx,:)];
     
    for j=1:size(countries,1)
        countries(j,3)=total(find(ismember(total(:,1),countries(j,1))==1),3);
        countries(j,4)=total(find(ismember(total(:,1),countries(j,1))==1),5);
        countries(j,5)=GDP_change_rate(find(ismember(GDP_change_rate(:,1),countries(j,1))==1),2);
    end
    % ===Peer filtering by performance===
         country_own=[total(find(ismember(total(:,1),country)==1),1),total(find(ismember(total(:,1),country)==1),3),total(find(ismember(total(:,1),country)==1),5),GDP_change_rate(find(ismember(GDP_change_rate(:,1),country)==1),2)];
         countries=sortrows(countries, 4);
         countries = countries(cell2mat(countries(:, 5)) > 0, :);
         countries = countries(cell2mat(countries(:, 4)) < cell2mat(country_own(:,3)), :);
         countries = countries(cell2mat(countries(:, 4)) < 0, :);
         countries(:,6)=num2cell(cell2mat(countries(:,5))./abs(cell2mat(countries(:,4))));
         countries=sortrows(countries, 6, 'descend');
     if isempty(countries)
         countries=[country_own(:,1),{0},country_own(:,2:end)];
     end

      average_intensity=repmat(cell2mat(countries(:,4)),1,5);
      average_GDP=repmat(cell2mat(countries(:,5)),1,5);
      average_gdppc=repmat(cell2mat(countries(:,2)),1,5);

      intensity= num2cell(cell2mat(intensity_final_2030(i,end)).* cumprod(1 +average_intensity,2));
      combine=[repmat(intensity_final_2030(i,:) ,size(intensity,1),1),intensity];
      intensity_final_2035=[intensity_final_2035;combine];
      gdppc=num2cell(cell2mat(gdppc_final_2030(i,end)).* cumprod(1 +average_gdppc,2));
      combine_gdppc=[repmat(gdppc_final_2030(i,:) ,size(gdppc,1),1),gdppc];
      gdppc_final_2035=[gdppc_final_2035;combine_gdppc];
    
      GDP=num2cell(cell2mat(GDP_final_2030(i,end)).* cumprod(1 +average_GDP,2));
      combine_GDP=[repmat(GDP_final_2030(i,:) ,size(GDP,1),1),GDP];
      GDP_final_2035=[GDP_final_2035;combine_GDP];
      CO2=num2cell(cell2mat(GDP_final_2035((size(GDP_final_2035,1)-size(countries,1)+1):end,23:27)).*cell2mat(intensity_final_2035((size(GDP_final_2035,1)-size(countries,1)+1):end,23:27)));
      combine_CO2=[repmat(CO2_final_2030(i,:) ,size(CO2,1),1),CO2];
      CO2_final_2035=[CO2_final_2035;combine_CO2];

end
intensity_final_2040=[];
gdppc_final_2040=[];
GDP_final_2040=[];
CO2_final_2040=[];
for i=1:size(gdppc_final_2035,1)
    country=gdppc_final_2035(i,1);
    target_2020=gdppc_final_2035(i,end);
    values = cell2mat(sorted_data(:,2)); 
    % ===Peer Selection===
    [~, closest_idx] = min(abs(values - cell2mat(target_2020))); 
    if closest_idx <= 5
        start_idx=1;
        end_idx=11;
    elseif closest_idx > total_rows - 5   
        start_idx=total_rows-10;
        end_idx=total_rows;
    else
    start_idx = closest_idx - 5; 
    end_idx = closest_idx + 5; 
    end
     countries=[sorted_data(start_idx:closest_idx-1,:);sorted_data(closest_idx+1:end_idx,:)];
    
    for j=1:size(countries,1)
        countries(j,3)=total(find(ismember(total(:,1),countries(j,1))==1),3);
        countries(j,4)=total(find(ismember(total(:,1),countries(j,1))==1),5);
        countries(j,5)=GDP_change_rate(find(ismember(GDP_change_rate(:,1),countries(j,1))==1),2);
    end
    % ===Peer filtering by performance===
         country_own=[total(find(ismember(total(:,1),country)==1),1),total(find(ismember(total(:,1),country)==1),3),total(find(ismember(total(:,1),country)==1),5),GDP_change_rate(find(ismember(GDP_change_rate(:,1),country)==1),2)];
         countries=sortrows(countries, 4);
         countries = countries(cell2mat(countries(:, 5)) > 0, :);
         countries = countries(cell2mat(countries(:, 4)) < cell2mat(country_own(:,3)), :);
         countries = countries(cell2mat(countries(:, 4)) < 0, :);
         countries(:,6)=num2cell(cell2mat(countries(:,5))./abs(cell2mat(countries(:,4))));
         countries=sortrows(countries, 6, 'descend');
     if isempty(countries)
         countries=[country_own(:,1),{0},country_own(:,2:end)];
     end

      average_intensity=repmat(cell2mat(countries(:,4)),1,5);
      average_GDP=repmat(cell2mat(countries(:,5)),1,5);
      average_gdppc=repmat(cell2mat(countries(:,2)),1,5);

      intensity= num2cell(cell2mat(intensity_final_2035(i,end)).* cumprod(1 +average_intensity,2));
      combine=[repmat(intensity_final_2035(i,:) ,size(intensity,1),1),intensity];
      intensity_final_2040=[intensity_final_2040;combine];
      gdppc=num2cell(cell2mat(gdppc_final_2035(i,end)).* cumprod(1 +average_gdppc,2));
      combine_gdppc=[repmat(gdppc_final_2035(i,:) ,size(gdppc,1),1),gdppc];
      gdppc_final_2040=[gdppc_final_2040;combine_gdppc];
    
      GDP=num2cell(cell2mat(GDP_final_2035(i,end)).* cumprod(1 +average_GDP,2));
      combine_GDP=[repmat(GDP_final_2035(i,:) ,size(GDP,1),1),GDP];
      GDP_final_2040=[GDP_final_2040;combine_GDP];
      CO2=num2cell(cell2mat(GDP_final_2040((size(GDP_final_2040,1)-size(countries,1)+1):end,28:32)).*cell2mat(intensity_final_2040((size(GDP_final_2040,1)-size(countries,1)+1):end,28:32)));
      combine_CO2=[repmat(CO2_final_2035(i,:) ,size(CO2,1),1),CO2];
      CO2_final_2040=[CO2_final_2040;combine_CO2];

end




xlswrite('F:\research\code_and_data\data_generation\excellence_scenarios.xlsx',CO2_final_2040,'all scenarios');
clearvars -except CO2_final_2040

country_names = CO2_final_2040(:, 1);

% Get unique country names
unique_countries = unique(country_names);

% Find max and min (range)
min_min = [];
max_max = [];

for i = 1:length(unique_countries)
    this_country = unique_countries{i};

    % Find all rows corresponding to the country
    rows = strcmp(country_names, this_country);

    
    numeric_data = cell2mat(CO2_final_2040(rows, 3:end));

    % Extract ISO code
    iso_code = CO2_final_2040{find(rows, 1), 2};

    % Compute column-wise minimum and maximum values
    
    min_col = min(numeric_data, [], 1);  % minimum of each column 
    max_col = max(numeric_data, [], 1);  % maximum of each column 

    % Save results (combine country name, ISO, and values)
    combine_1 = [this_country, iso_code, num2cell(max_col)];
    combine_2 = [this_country, iso_code, num2cell(min_col)];
   

    % Append to final results
    min_min = [min_min; combine_2];
    max_max = [max_max; combine_1];
end



xlswrite('F:\research\code_and_data\data_generation\excellence_scenarios.xlsx',max_max,'max');
xlswrite('F:\research\code_and_data\data_generation\excellence_scenarios.xlsx',min_min,'min');

