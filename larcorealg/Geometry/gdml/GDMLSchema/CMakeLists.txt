# install schema files
file(GLOB gdml_schema *.xsd)
install( FILES ${gdml_schema} DESTINATION ${${PROJECT_NAME}_GDML_DIR}/GDMLSchema )
