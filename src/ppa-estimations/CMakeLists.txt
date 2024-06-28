# ---------------------------------------------------------------------------- #
# source 
# ---------------------------------------------------------------------------- #
target_sources(${ILANG_LIB_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/ppa_main.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ppa_hardware_block.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ppa_model_registrar.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ppa_register_models.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ppa_callbacks.cc

)

