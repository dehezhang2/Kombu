#configuration
CONFIG +=  shared def_files_disabled exceptions no_mocdepend release stl qt_framework  x86_64
QT_ARCH = macosx
QT_EDITION = OpenSource
QT_CONFIG +=  minimal-config small-config medium-config large-config full-config no-pkg-config qt3support accessibility opengl shared reduce_exports ipv6 getaddrinfo ipv6ifname getifaddrs png no-freetype system-zlib nis cups iconv openssl-linked corewlan concurrent xmlpatterns multimedia audio-backend svg script scripttools declarative release qt_framework  x86_64

#versioning
QT_VERSION = 4.8.7
QT_MAJOR_VERSION = 4
QT_MINOR_VERSION = 8
QT_PATCH_VERSION = 7

#namespaces
QT_LIBINFIX = 
QT_NAMESPACE = 
QT_NAMESPACE_MAC_CRC = 

QMAKE_RPATHDIR += "/usr/local/Cellar/qt@4/4.8.7_6/lib"
