esma_set_this ()

set (SRCS
  ACCoeff_Binary_IO.f90
  ACCoeff_Define.f90
  ADA_Module.f90
  AOvar_Define.f90
  ASvar_Define.f90
  AerosolCoeff_Binary_IO.f90
  AerosolCoeff_Define.f90
  AerosolCoeff_netCDF_IO.f90
  Azimuth_Emissivity_F6_Module.f90
  Azimuth_Emissivity_Module.f90
  Binary_File_Utility.f90
  CRTM_AOD_Module.f90
  CRTM_Adjoint_Module.f90
  CRTM_AerosolCoeff.f90
  CRTM_AerosolScatter.f90
  CRTM_Aerosol_Define.f90
  CRTM_AncillaryInput_Define.f90
  CRTM_AntennaCorrection.f90
  CRTM_AtmAbsorption.f90
  CRTM_AtmOptics.f90
  CRTM_AtmOptics_Define.f90
  CRTM_Atmosphere.f90
  CRTM_Atmosphere_Define.f90
  CRTM_ChannelInfo_Define.f90
  CRTM_CloudCoeff.f90
  CRTM_CloudScatter.f90
  CRTM_Cloud_Define.f90
  CRTM_CloudCover_Define.f90
  CRTM_Fastem1.f90
  CRTM_FastemX.f90
  CRTM_Forward_Module.f90
  CRTM_GeometryInfo.f90
  CRTM_GeometryInfo_Define.f90
  CRTM_Geometry_Define.f90
  CRTM_IRSSEM.f90
  CRTM_IR_Ice_SfcOptics.f90
  CRTM_IR_Land_SfcOptics.f90
  CRTM_IR_Snow_SfcOptics.f90
  CRTM_IR_Water_SfcOptics.f90
  CRTM_IRiceCoeff.f90
  CRTM_IRlandCoeff.f90
  CRTM_IRsnowCoeff.f90
  CRTM_IRwaterCoeff.f90
  CRTM_Interpolation.f90
  CRTM_K_Matrix_Module.f90
  CRTM_LifeCycle.f90
  CRTM_LowFrequency_MWSSEM.f90
  CRTM_MW_Ice_SfcOptics.f90
  CRTM_MW_Land_SfcOptics.f90
  CRTM_MW_Snow_SfcOptics.f90
  CRTM_MW_Water_SfcOptics.f90
  CRTM_MWwaterCoeff.f90
  CRTM_Model_Profiles.f90
  CRTM_MoleculeScatter.f90
  CRTM_NLTECorrection.f90
  CRTM_Options_Define.f90
  CRTM_Parameters.f90
  CRTM_Planck_Functions.f90
  CRTM_Predictor.f90
  CRTM_Predictor_Define.f90
  CRTM_RTSolution.f90
  CRTM_RTSolution_Define.f90
  CRTM_SEcategory.f90
  CRTM_SensorData_Define.f90
  CRTM_SensorInfo.f90
  CRTM_SfcOptics.f90
  CRTM_SfcOptics_Define.f90
  CRTM_SpcCoeff.f90
  CRTM_Surface_Define.f90
  CRTM_Tangent_Linear_Module.f90
  CRTM_TauCoeff.f90
  CRTM_Utility.f90
  CRTM_VIS_Ice_SfcOptics.f90
  CRTM_VIS_Land_SfcOptics.f90
  CRTM_VIS_Snow_SfcOptics.f90
  CRTM_VIS_Water_SfcOptics.f90
  CRTM_VISiceCoeff.f90
  CRTM_VISlandCoeff.f90
  CRTM_VISsnowCoeff.f90
  CRTM_VISwaterCoeff.f90
  CSvar_Define.f90
  CloudCoeff_Binary_IO.f90
  CloudCoeff_Define.f90
	CloudCoeff_netCDF_IO.f90
  Common_RTSolution.f90
  Compare_Float_Numbers.f90
  DateTime_Utility.f90
  Date_Utility.f90
  Ellison.f90
  Emission_Module.f90
  Endian_Utility.f90
  File_Utility.f90
  FitCoeff_Define.f90
  Foam_Utility_Module.f90
  Fresnel.f90
  Fundamental_Constants.f90
  Guillou.f90
  Hyperbolic_Step.f90
  IRwaterCoeff_Define.f90
  LSEatlas_Define.f90
  Large_Scale_Correction_Module.f90
  Liu.f90
  MWwaterCoeff_Define.f90
  MWwaterLUT_Define.f90
  Message_Handler.f90
  NESDIS_AMSRE_SICEEM_Module.f90
  NESDIS_AMSRE_SNOWEM_Module.f90
  NESDIS_AMSU_SICEEM_Module.f90
  NESDIS_AMSU_SnowEM_Module.f90
  NESDIS_ATMS_SnowEM_Module.f90
  NESDIS_ATMS_SeaICE_LIB.f90
  NESDIS_ATMS_SeaICE_Module.f90
  NESDIS_LandEM_Module.f90
  NESDIS_MHS_SICEEM_Module.f90
  NESDIS_MHS_SnowEM_Module.f90
  NESDIS_SEAICE_PHYEM_MODULE.f90
  NESDIS_SSMIS_SeaIceEM_Module.f90
  NESDIS_SSMIS_SnowEM_Module.f90
  NESDIS_SSMI_Module.f90
  NESDIS_SSMI_SIceEM_Module.f90
  NESDIS_SSMI_SnowEM_Module.f90
  NESDIS_SnowEM_ATMS_Parameters.f90
  NESDIS_SnowEM_Parameters.f90
  NLTECoeff_Binary_IO.f90
  NLTECoeff_Define.f90
  NLTE_Parameters.f90
  NLTE_Predictor_Define.f90
  NLTE_Predictor_IO.f90
  ODAS_AtmAbsorption.f90
  ODAS_Binary_IO.f90
  ODAS_Define.f90
  ODAS_Predictor.f90
  ODAS_Predictor_Define.f90
  ODAS_TauCoeff.f90
  ODPS_AtmAbsorption.f90
  ODPS_Binary_IO.f90
  ODPS_CoordinateMapping.f90
  ODPS_Define.f90
  ODPS_Predictor.f90
  ODPS_Predictor_Define.f90
  ODPS_TauCoeff.f90
  ODSSU_AtmAbsorption.f90
  ODSSU_Binary_IO.f90
  ODSSU_Define.f90
  ODSSU_TauCoeff.f90
  ODZeeman_AtmAbsorption.f90
  ODZeeman_Predictor.f90
  ODZeeman_TauCoeff.f90
  PAFV_Define.f90
  Profile_Utility_Parameters.f90
  RTV_Define.f90
  Reflection_Correction_Module.f90
  SEcategory_Define.f90
  SOI_Module.f90
  SSU_Input_Define.f90
  Search_Utility.f90
  SensorInfo_Parameters.f90
  Slope_Variance.f90
  Small_Scale_Correction_Module.f90
  Sort_Utility.f90
  SpcCoeff_Binary_IO.f90
  SpcCoeff_Define.f90
  Spectral_Units_Conversion.f90
  String_Utility.f90
  Subset_Define.f90
  TauCoeff_Define.f90
  Timing_Utility.f90
  Type_Kinds.f90
  UnitTest_Define.f90
  Zeeman_Input_Define.f90
  Zeeman_Utility.f90
  iAtm_Define.f90
  CRTM_Module.F90
)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set (CMAKE_Fortran_FLAGS_RELEASE "-free ${FOPT2} ${LITTLE_ENDIAN} ${BYTERECLEN}")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(WARNING "Intel Compiler and our debugging flags have issues with ADA_Module.f90. So for now we turn off checking compiling ADA_Module.f90")
    set_source_files_properties(ADA_Module.f90 PROPERTIES COMPILE_OPTIONS -nocheck)
  endif ()
endif ()


esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES NetCDF::NetCDF_Fortran
  )

set_target_properties(${this} PROPERTIES POSITION_INDEPENDENT_CODE OFF)
