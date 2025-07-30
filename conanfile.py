from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.files import copy
import os

required_conan_version = ">=2.0.9"

class SvgfillConan(ConanFile):
    name = "svgfill"
    version = "0.1.0"
    description = "An application to fill areas bounded by unconnected lines in SVG"
    license = "LGPL-2.1-or-later"
    url = "https://github.com/conan-io/conan-center-index"
    homepage = "https://github.com/IfcOpenShell/svgfill"
    topics = ("svg", "ifc")
    settings = "os", "arch", "compiler", "build_type"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
    }
    default_options = {
        "shared": False,
        "fPIC": True,
    }
    implements = ["auto_shared_fpic"]

    exports_sources = "CMakeLists.txt", "src/**"
    
    # no exports_sources attribute, but export_sources(self) method instead
    # def export_sources(self):
    #     export_conandata_patches(self)

    def layout(self):
        cmake_layout(self, src_folder="")

    def requirements(self):
        self.requires("svgpp/1.3.1")
        self.requires("cgal/6.0.1")
        self.requires("boost/1.88.0", override=True, transitive_headers=True)
        self.requires("libxml2/2.13.8")
        self.requires("nlohmann_json/3.12.0")
        self.requires("gmp/6.3.0")
        self.requires("mpfr/4.2.1")

    def validate(self):
        check_min_cppstd(self, 17)

    def build_requirements(self):
        self.tool_requires("cmake/[>=3.16 <4]")

    # def source(self):
    #     get(self, **self.conan_data["sources"][self.version], strip_root=True)
    #     # Using patches is always the last resort to fix issues. If possible, try to fix the issue in the upstream project.
    #     apply_conandata_patches(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["svgfill"]
        