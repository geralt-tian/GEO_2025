add_library(LIB_BERT 
    he.cpp
    nonlinear.cpp
    linear.cpp
    bert_utils.cpp
    bert.cpp
    )

target_link_libraries(LIB_BERT SCI-FloatingPoint)
target_link_libraries(LIB_BERT SCI-HE)

# add_executable(BERT_NONLINEAR "test_nonlinear.cpp")
# target_link_libraries(BERT_NONLINEAR NONLINEAR)

add_executable(BOLT_BERT "bolt_bert.cpp")
target_link_libraries(BOLT_BERT LIB_BERT)


# add_executable(BOLT_WE_BERT "bolt_bert_word_elimination.cpp")
# target_link_libraries(BOLT_WE_BERT LIB_BERT)

# OUR SOFTMAX VERSION - Using softmax_fix_our
add_library(LIB_BERT_OUR
    he.cpp
    nonlinear_our.cpp
    linear.cpp
    bert_utils.cpp
    bert_our.cpp
    )

target_link_libraries(LIB_BERT_OUR SCI-FloatingPoint)
target_link_libraries(LIB_BERT_OUR SCI-HE)

add_executable(BOLT_BERT_OUR "bolt_bert_our.cpp")
target_link_libraries(BOLT_BERT_OUR LIB_BERT_OUR)
