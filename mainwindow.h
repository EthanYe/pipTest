#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include<QPixmap>
#include <QMainWindow>
#include<QLabel>
#include<QActionGroup>
#include<QButtonGroup>
#include<QPushButton>
#include<QPainter>
#include<pipsys.h>
#include<QMetaEnum>
#include<QMap>
#include<QRadioButton>
QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    Pipsys *ps;
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    //action
     QAction *openAction;

    QButtonGroup *reservoirBtns;
    QButtonGroup *pipeBtns;
    QButtonGroup *nodeBtns;
    QButtonGroup *branchBtns;
    QButtonGroup *elementsBtns;
QPushButton* virtualBtn;
    QActionGroup *elementActs;
    int nRe=0;
    void mousePressEvent(QMouseEvent *event);
    QMetaEnum elementEnum;
    QStringList elementsList;
    QList<QButtonGroup *> eleBGList;
    void addElementAction();
    QStringList button_styleList;
    QList<int > eleCount;


    void reservoirDialog(int id);
};
#endif // MAINWINDOW_H
