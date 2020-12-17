QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    QPipe.cpp \
    elements.cpp \
    hydroment.cpp \
    main.cpp \
    mainwindow.cpp \
    myQT.cpp \
    my_math.cpp \
    pipsys.cpp

HEADERS += \
    Qpipe.h \
    elements.h \
    hydroment.h \
    mainwindow.h \
    myQT.h \
    my_math.h \
    pipsys.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    res.qrc
INCLUDEPATH+= E:/DeskFile/2018C/third_party/armadillo/include\

LIBS += \
E:\DeskFile\2018C\third_party\armadillo\lib_win64\libopenblas.lib\
