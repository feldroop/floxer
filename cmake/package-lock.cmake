CPMDeclarePackage (
    sharg
    NAME sharg
    GITHUB_REPOSITORY seqan/sharg-parser
    GIT_TAG 39f65a4890f8c5108af2b5c7974893ff6ed87e50
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
    OPTIONS
        "SHARG_NO_TDL ON"
)

CPMDeclarePackage (
    spdlog
    NAME spdlog
    GITHUB_REPOSITORY gabime/spdlog
    GIT_TAG v1.13.0
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
)

CPMDeclarePackage (
    IVio
    NAME IVio
    GITHUB_REPOSITORY iv-project/IVio
    GIT_TAG 8fcb901bc08dd48d050bcfbcf73d5fe40db6c90a
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
)

CPMDeclarePackage (
    IVSigma
    NAME IVSigma
    GITHUB_REPOSITORY iv-project/IVSigma
    GIT_TAG 4b76036ed331b91b7a08fc7249ac51b4def71c74
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
)

CPMDeclarePackage (
    fmindex-collection
    NAME fmindex-collection
    GITHUB_REPOSITORY SGSSGene/fmindex-collection
    GIT_TAG b0e311f122f6ae9d15d470ca88ac21f028dc010a
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
)

CPMDeclarePackage (
    cereal
    NAME cereal
    GITHUB_REPOSITORY USCiLab/cereal
    GIT_TAG v1.3.2
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
    OPTIONS
        "SKIP_PERFORMANCE_COMPARISON ON"
        "BUILD_SANDBOX OFF"
        "BUILD_DOC OFF"
)

CPMDeclarePackage (
    seqan3
    NAME seqan3
    GITHUB_REPOSITORY seqan/seqan3
    GIT_TAG 41a17ad84058e3c771efe0f2d10764c31f121583
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
    OPTIONS
        "INSTALL_SEQAN3 OFF"
        "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

CPMDeclarePackage (
    interval-tree
    NAME interval-tree
    GITHUB_REPOSITORY 5cript/interval-tree
    GIT_TAG 4b91cf434eebceef4b8509c97b568d7efef14b20
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
)

CPMDeclarePackage (
    BS_thread_pool
    NAME BS_thread_pool
    GITHUB_REPOSITORY bshoshany/thread-pool
    VERSION 4.1.0
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
)

CPMDeclarePackage (
    cpp-channel
    NAME cpp-channel
    GITHUB_REPOSITORY andreiavrammsd/cpp-channel
    GIT_TAG 38ffdec0b7eec2acbbe98026663d1a7e1acc29eb
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
)

CPMDeclarePackage (
    use_ccache
    NAME use_ccache
    GITHUB_REPOSITORY seqan/cmake-scripts
    GIT_TAG d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
    SOURCE_SUBDIR ccache
)

CPMDeclarePackage(
    googletest
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    GIT_TAG 5197b1a8e6a1ef9f214f4aa537b0be17cbf91946
    EXCLUDE_FROM_ALL TRUE
    SYSTEM TRUE
    OPTIONS
        "INSTALL_GTEST OFF"
        "gtest_force_shared_crt"
)