load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "com_google_googletest",
    strip_prefix = "googletest-5ab508a01f9eb089207ee87fd547d290da39d015",
    urls = ["https://github.com/google/googletest/archive/5ab508a01f9eb089207ee87fd547d290da39d015.zip"],
)

git_repository(
    name = "happly",
    branch = "master",
    build_file_content = """
cc_library(
    name = "happly",
    srcs = [],
    includes = ['.'],
    hdrs = ["happly.h"],
    visibility = ["//visibility:public"],
)
""",
    remote = "https://github.com/nmwsharp/happly.git",
)
