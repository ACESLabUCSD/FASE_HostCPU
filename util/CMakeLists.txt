add_library (Util_Log OBJECT
  "log.cpp")

add_library (Util_Util OBJECT
  "util.cpp")

add_executable(Util_Test
  "util_test.cpp"
  "$<TARGET_OBJECTS:Util_Util>"
  "$<TARGET_OBJECTS:Util_Log>")
if (Boost_FOUND)
	target_link_libraries(Util_Test ${Boost_LIBRARIES})
endif (Boost_FOUND)
if (OPENSSL_FOUND)
	target_link_libraries(Util_Test ${OPENSSL_LIBRARIES})
endif (OPENSSL_FOUND)
#target_link_libraries(Util_Test -lrt)

add_executable(txt2bin
  "txt2bin.cpp"
  "$<TARGET_OBJECTS:Garbled_Circuit_Garbled_Circuit>"
  "$<TARGET_OBJECTS:Garbled_Circuit_Garbled_Circuit_Low_Mem>"
  "$<TARGET_OBJECTS:Crypto_BN>"
  "$<TARGET_OBJECTS:Crypto_OT>"
  "$<TARGET_OBJECTS:Crypto_OT_Extension>"
  "$<TARGET_OBJECTS:SCD_SCD>"
  "$<TARGET_OBJECTS:Util_Util>"
  "$<TARGET_OBJECTS:TCPIP_TCPIP>"
  "$<TARGET_OBJECTS:Util_Log>")

if (OPENSSL_FOUND)
	target_link_libraries(txt2bin ${OPENSSL_LIBRARIES})
endif (OPENSSL_FOUND)
if (Boost_FOUND)
  target_link_libraries(txt2bin ${Boost_LIBRARIES})
endif (Boost_FOUND)
target_link_libraries(txt2bin -lrt)

add_test (Util_Util_Test
  Util_Test --log2std)

set_tests_properties (Util_Util_Test
	PROPERTIES
	FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;failed")
