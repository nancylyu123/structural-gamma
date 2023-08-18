
import PatientGamma_Functions as PG
import numpy as np
import pydicom
import pymedphys

##############################################
# Adjustable Parameters

# List of Patients/Structures for which to calculate Gamma Calcs
PatientDetails_filename = 'filename.csv'

# List of Stuctures without a Tolerance for which to calculate Gamma
NonToleranceStructureList = ['95% ISODOSE', 'PTV', 'PTV1', 'PTV2', 'PTV3', 'PTV4', 'TOTAL']

# Calaculate Gamma for structures with a tolerance
ProcessToleranceStructures = True

# Gamma Parameters
dd               = 2      # Dose difference percent
dta              = 2      # Distance to agreement
LowDoseThreshold = 0.0001 # Low dose cuttoff percent
InterpFrac       = 10     # Interpolation fraction
MaxGamma         = 2      # Maximum Value to calculate Gamma to

# Save Gamma Index arrays for later plotting
#   (note: if ProcessToleranceStructures is False then Gamma Index arrays 
#    will NOT be saved as they would be empty)
SaveDataArrays = True

##############################################

# Load the Standard Structure Tolerances
StructureTolerance_dictionary = PG.LoadStructureTolerance()

# Load the Patient Data into a list of Patient objects
PatientList = PG.LoadPatients(PatientDetails_filename)

# Open file for comma seperated results
f_results = open(PG.Results_filename, 'w')
# Write results headings
f_results.write("Patient,ID,Prescription Dose,Structure,Standard Structure,Tolerance,Max Body or External Dose,Max Structure Dose,")
f_results.write("Global,# Global Evaluated,")
f_results.write("Local,# Local Evaluated,")
f_results.write("QUANTEC Normalised,# QUANTEC Evaluated,")
f_results.write("RO Normalised,# RO Evaluated,RO Dose Tolerance,")
f_results.write("TPS Volume,MC Volume,")
f_results.write("TPS Mean,MC Mean,Mean Diff Global,Mean Diff Local,Mean Diff QUANTEC Normalised,Mean Diff RO Normalised,")
f_results.write("Ref file name,Eval file name,")
f_results.write("Dose Difference %,Distance to Agreement,Low Dose Threshold,Interpolation Fraction,Max Gamma Calculation")
f_results.write("\n")

for Patient in PatientList:
    print(f'Processing Patient {Patient.Name}')

    # Get Prescition Dose
    PrescriptionDose = 0
    RP_filename = f'{PG.DoseFile_dirname}/{Patient.Name}_{Patient.ID}/RP.{Patient.ID}.{Patient.Plan}.dcm'
    RP_dicom    = pydicom.dcmread(RP_filename)
    PrescriptionDose = float(RP_dicom.DoseReferenceSequence[0].DeliveryMaximumDose)

    # Get TPS Dose for either Body or External
    TPS_Dose_Body_filename       = f'{PG.DoseFile_dirname}/{Patient.Name}_{Patient.ID}/RD.{Patient.ID}.Dose_{Patient.Plan}.dcm_{Patient.Body}.dcm'
    Dicom_TPS_Dose_Body          = pydicom.dcmread(TPS_Dose_Body_filename,  force=True)
    TPS_Axes_Body, TPS_Dose_Body = pymedphys.dicom.zyx_and_dose_from_dataset(Dicom_TPS_Dose_Body)
    (z_TPS, y_TPS, x_TPS)        = TPS_Axes_Body
    MaxTPSBodyDose               = TPS_Dose_Body.max()
    
    # Get Evaluation dataset (MC Body or External)
    MC_Dose_Body_filename        = f'{PG.DoseFileMC_dirname}/{Patient.Name}_{Patient.ID}/RD.{Patient.ID}.Dose_{Patient.PlanMC}.dcm_{Patient.Body}.dcm'
    Dicom_MC_Dose_Body           = pydicom.dcmread(MC_Dose_Body_filename,  force=True)
    MC_Axes_Body,  MC_Dose_Body  = pymedphys.dicom.zyx_and_dose_from_dataset(Dicom_MC_Dose_Body)
    (z_MC,  y_MC,  x_MC)         = MC_Axes_Body
    
    # Check TPS and MC Dose files for an offset and resample the MC Dose file if necessary
    MC_Dose_Aligned = PG.DoseAlign(TPS_Axes_Body, TPS_Dose_Body, MC_Axes_Body, MC_Dose_Body)

    # Calculate voxel size in cc
    VoxelsPerCC = 1000 / ((z_TPS[1]-z_TPS[0]) * (y_TPS[1]-y_TPS[0]) * (x_TPS[1]-x_TPS[0]))

    # List of Gamma Index arrays created for each normalisation method
    Gamma_list_Global            = []
    Gamma_list_Local             = []
    Gamma_list_Tolerance         = []
    Gamma_list_ROTolerance       = []
    
    # List of all Patient structures processed
    StructureName_list           = []
    
    # Calculate Gamma for each Patient Structure
    for Structures in Patient.Structures:
        Structure     = Structures[0]
        PlanStructure = Structures[1]
        ROTolerance   = Structures[2]
        print(f"Processing structure {PlanStructure}\n")

        if len(PlanStructure) > 0:

            # Get the TPS/Reference dataset for the Structure
            TPS_Dose_filename     = f'{PG.DoseFile_dirname}/{Patient.Name}_{Patient.ID}/RD.{Patient.ID}.Dose_{Patient.Plan}.dcm_{PlanStructure}.dcm'
            Dicom_TPS_Dose        = pydicom.dcmread(TPS_Dose_filename, force=True)
            TPS_Axes, TPS_Dose    = pymedphys.dicom.zyx_and_dose_from_dataset(Dicom_TPS_Dose)
            (z_TPS, y_TPS, x_TPS) = TPS_Axes
            
            # Get the MC/Evaluation dataset for the Structure
            MC_Dose_filename      = f'{PG.DoseFileMC_dirname}/{Patient.Name}_{Patient.ID}/RD.{Patient.ID}.Dose_{Patient.PlanMC}.dcm_{PlanStructure}.dcm'
            Dicom_MC_Dose         = pydicom.dcmread(MC_Dose_filename, force=True)
            MC_Axes, MC_Dose      = pymedphys.dicom.zyx_and_dose_from_dataset(Dicom_MC_Dose)
            (z_MC,  y_MC,  x_MC)  = MC_Axes
            
            # Set the MC Structure as a mask for the TPS (i.e. set the TPS Dose to 0 where the MC Dose is 0)
            TPS_Dose[np.where(MC_Dose_Aligned==0)] = 0

            MaxTPSStructureDose = TPS_Dose.max()
                
            # Determine if this is a Standard Structure
            if Structure.upper() in StructureTolerance_dictionary.keys() and ProcessToleranceStructures:
                
                # Get QUANTEC dose tolerance
                DoseTolerance = StructureTolerance_dictionary[Structure.upper()]

                # Write values to results file
                f_results.write(f"{Patient.Name},{Patient.ID},{PrescriptionDose},{Structure},{PlanStructure},{DoseTolerance},{MaxTPSBodyDose},{MaxTPSStructureDose},")
                
                ### Perform Gamma Analysis
                
                ### Global 2% 2mm
                Gamma, Gamma_calculated, Gamma_pass_ratio = PG.CalculateGamma(False, MaxTPSBodyDose, TPS_Axes, TPS_Dose, MC_Axes_Body, MC_Dose_Body, dd, dta, LowDoseThreshold, InterpFrac, MaxGamma)
                f_results.write(str("{0:.2f}".format(Gamma_pass_ratio*100)) + "," + str(len(Gamma_calculated)) + ",")
                Gamma_list_Global.append(Gamma)

                ### Local 2% 2mm
                Gamma, Gamma_calculated, Gamma_pass_ratio = PG.CalculateGamma(True, None, TPS_Axes, TPS_Dose, MC_Axes_Body, MC_Dose_Body, dd, dta, LowDoseThreshold, InterpFrac, MaxGamma)
                f_results.write(str("{0:.2f}".format(Gamma_pass_ratio*100)) + "," + str(len(Gamma_calculated)) + ",")
                Gamma_list_Local.append(Gamma)
                
                ### Tolerance Normalised 2% 2mm
                if DoseTolerance > 0:
                    Gamma, Gamma_calculated, Gamma_pass_ratio = PG.CalculateGamma(False, DoseTolerance, TPS_Axes, TPS_Dose, MC_Axes_Body, MC_Dose_Body, dd, dta, LowDoseThreshold, InterpFrac, MaxGamma)
                    f_results.write(str("{0:.2f}".format(Gamma_pass_ratio*100)) + "," + str(len(Gamma_calculated)) + ",")
                    Gamma_list_Tolerance.append(Gamma)
                else:
                    f_results.write(",,")
                
                ### RO Tolerance Normalised 2% 2mm
                if ROTolerance > 0:
                    Gamma, Gamma_calculated, Gamma_pass_ratio = PG.CalculateGamma(False, ROTolerance, TPS_Axes, TPS_Dose, MC_Axes_Body, MC_Dose_Body, dd, dta, LowDoseThreshold, InterpFrac, MaxGamma)
                    f_results.write(str("{0:.2f}".format(Gamma_pass_ratio*100)) + "," + str(len(Gamma_calculated)) + "," + str("{0:.1f}".format(ROTolerance)) + ",")
                    Gamma_list_ROTolerance.append(Gamma)
                else:
                    f_results.write(",,,")
                      
                
                StructureName_list.append(Structure)
                
                
                ### Calculate Structure Volume Statistics
                
                # Sort TPS Dose values
                TPS_Dose_flat = np.sort(TPS_Dose[TPS_Dose>0], axis=None)
                DataPoints    = TPS_Dose_flat.size
                TPS_RelativeVolume = np.linspace(DataPoints, 0, num=DataPoints) / DataPoints

                # Sort MC Dose values
                MC_Dose_flat = np.sort(MC_Dose[MC_Dose>0], axis=None)
                DataPoints = MC_Dose_flat.size
                MC_RelativeVolume = np.linspace(DataPoints, 0, num=DataPoints) / DataPoints

                # Write Volume
                f_results.write(str("{0:.4f}".format(TPS_Dose_flat.size/VoxelsPerCC)) + ",")
                f_results.write(str("{0:.4f}".format(MC_Dose_flat.size/VoxelsPerCC)) + ",")

                # Calculate Mean
                TPS_Mean = TPS_Dose_flat.mean()
                MC_Mean  = MC_Dose_flat.mean()
                
                # Write Mean values
                f_results.write(str("{0:.4f}".format(TPS_Mean)) + "," + str("{0:.4f}".format(MC_Mean)) + ",")
                                
                # Write Mean Difference (Globaly Normalised)
                MeanDiff = (MC_Mean - TPS_Mean) / MaxTPSBodyDose * 100
                f_results.write(str("{0:.4f}".format(MeanDiff)) + ",")
    
                # Write Mean Difference (Localy Normalised)
                MeanDiff = (MC_Mean - TPS_Mean) / TPS_Mean * 100
                f_results.write(str("{0:.4f}".format(MeanDiff)) + ",")
                
                if DoseTolerance > 0:
                    # Write Mean Difference (QUANTEC Normalised)
                    MeanDiff = (MC_Mean - TPS_Mean) / DoseTolerance * 100
                    f_results.write(str("{0:.4f}".format(MeanDiff)) + ",")
                else:
                    f_results.write(",")
                    
                if ROTolerance > 0:
                    # Write Mean Difference (RO Normalised)
                    MeanDiff = (MC_Mean - TPS_Mean) / ROTolerance * 100
                    f_results.write(str("{0:.4f}".format(MeanDiff)) + ",")
                else:
                    f_results.write(",")
                    
                # Add file names and Gamma parameters to the output
                f_results.write(f"{TPS_Dose_filename},{MC_Dose_Body_filename},")
                f_results.write(f"{dd},{dta},{LowDoseThreshold},{InterpFrac},{MaxGamma}\n")
                f_results.flush()

            else:
                # Perform Gama Analysis for Structures without a Tolerance
    
                if (Structure.upper() in NonToleranceStructureList):
    
                    ### Perform Gamma Analysis for
                    ### Global              2% 2mm
                    ### Local               2% 2mm
    
                    f_results.write(f"{Patient.Name},{Patient.ID},{PrescriptionDose},{Structure},{PlanStructure},0,{MaxTPSBodyDose},{MaxTPSStructureDose},")
                
                    ### Global 2% 2mm TH5
                    Gamma, Gamma_calculated, Gamma_pass_ratio = PG.CalculateGamma(False, MaxTPSBodyDose, TPS_Axes, TPS_Dose, MC_Axes_Body, MC_Dose_Body, dd, dta, LowDoseThreshold, InterpFrac, MaxGamma)
                    f_results.write(str("{0:.2f}".format(Gamma_pass_ratio*100)) + "," + str(len(Gamma_calculated)) + ",")
    
                    if Structure.upper() == 'TOTAL':
                        TotalGlobalGamma = Gamma

                    ### Local 2% 2mm TH5
                    Gamma, Gamma_calculated, Gamma_pass_ratio = PG.CalculateGamma(True, None, TPS_Axes, TPS_Dose, MC_Axes_Body, MC_Dose_Body, dd, dta, LowDoseThreshold, InterpFrac, MaxGamma)
                    f_results.write(str("{0:.2f}".format(Gamma_pass_ratio*100)) + "," + str(len(Gamma_calculated)) + ",")

                    if Structure.upper() == 'TOTAL':
                        TotalLocalGamma = Gamma
                        
                    f_results.write(",,,,,")


                    ### Calculate Structure Volume Statistics
                    
                    # Sort TPS Dose values
                    TPS_Dose_flat = np.sort(TPS_Dose[TPS_Dose>0], axis=None)
                    DataPoints    = TPS_Dose_flat.size
                    TPS_RelativeVolume = np.linspace(DataPoints, 0, num=DataPoints) / DataPoints
    
                    # Sort MC Dose values
                    MC_Dose_flat = np.sort(MC_Dose[MC_Dose>0], axis=None)
                    DataPoints = MC_Dose_flat.size
                    MC_RelativeVolume = np.linspace(DataPoints, 0, num=DataPoints) / DataPoints
    
                    # Write Volume assuming 2mm Voxels 1cc is 125 voxels
                    f_results.write(str("{0:.4f}".format(TPS_Dose_flat.size/VoxelsPerCC)) + ",")
                    f_results.write(str("{0:.4f}".format(MC_Dose_flat.size/VoxelsPerCC)) + ",")

                    # Calculate Mean values
                    TPS_Mean = TPS_Dose_flat.mean()
                    MC_Mean  = MC_Dose_flat.mean()
                    
                    # Write Mean values
                    f_results.write(str("{0:.4f}".format(TPS_Mean)) + "," + str("{0:.4f}".format(MC_Mean)) + ",")
                                    
                    # Write Mean Difference (Globaly Normalised)
                    MeanDiff = (MC_Mean - TPS_Mean) / MaxTPSBodyDose * 100
                    f_results.write(str("{0:.4f}".format(MeanDiff)) + ",")
        
                    # Write Mean Difference (Localy Normalised)
                    MeanDiff = (MC_Mean - TPS_Mean) / TPS_Mean * 100
                    f_results.write(str("{0:.4f}".format(MeanDiff)) + ",")
                    
                    f_results.write(",,")
    
                    # Add file names and Gamma parameters to the output
                    f_results.write(f"{TPS_Dose_filename},{MC_Dose_Body_filename},")
                    f_results.write(f"{dd},{dta},{LowDoseThreshold},{InterpFrac},{MaxGamma}\n")
                    f_results.flush()
                    
                else:
                    f_results.write(f"{Patient.Name},{Patient.ID},{PrescriptionDose},{Structure},{PlanStructure},\n")
                    f_results.flush()

        else:
            f_results.write(f"{Patient.Name},{Patient.ID},{PrescriptionDose},{Structure},\n")
            f_results.flush()

    f_results.write("\n")
    
         
    #######################################################################
    ## Save numpy arrays for later processing
    
    # if SaveDataArrays and ProcessToleranceStructures:
    if SaveDataArrays:
        print(f'Combining Gamma Arrays for {Patient.Name}')
        
        def CombineGamma(Gamma_List):
            # Combine Gamma Index arrays which contain nan values
            CombinedGamma    = np.empty_like(TPS_Dose_Body)
            CombinedGamma[:] = np.nan
            for Gamma in Gamma_List:
                CombinedGamma = np.where(np.isnan(CombinedGamma), Gamma, CombinedGamma)
            return CombinedGamma
        
        # Get combined numpy Gamma Index arrays and Hot/Cold arrays
        Gamma_Global    = CombineGamma(Gamma_list_Global)
        Gamma_Local     = CombineGamma(Gamma_list_Local)
        Gamma_Tolerance = CombineGamma(Gamma_list_Tolerance)
        if Gamma_list_ROTolerance:
            Gamma_ROTolerance = CombineGamma(Gamma_list_ROTolerance)

        print(f'Saving Gamma Arrays for {Patient.Name}')
        
        # Save the numpy arrays
        np.save(f'{PG.GammaArrays_dirname}/{Patient.Name}_{Patient.ID}_Global.npy', Gamma_Global)
        np.save(f'{PG.GammaArrays_dirname}/{Patient.Name}_{Patient.ID}_Local.npy', Gamma_Local)
        np.save(f'{PG.GammaArrays_dirname}/{Patient.Name}_{Patient.ID}_Tolerance.npy', Gamma_Tolerance)
        if Gamma_list_ROTolerance:
            np.save(f'{PG.GammaArrays_dirname}/{Patient.Name}_{Patient.ID}_ROTolerance.npy', Gamma_ROTolerance)

    # Save Total Gamma Index values
    if SaveDataArrays and 'TOTAL' in NonToleranceStructureList:
        # Save the numpy array for Total Gamma
        np.save(f'{PG.GammaArrays_dirname}/{Patient.Name}_{Patient.ID}_Total_Global.npy', TotalGlobalGamma)
        np.save(f'{PG.GammaArrays_dirname}/{Patient.Name}_{Patient.ID}_Total_Local.npy', TotalLocalGamma)

f_results.close()