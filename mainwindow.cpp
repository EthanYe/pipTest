#include "mainwindow.h"
#include "ui_mainwindow.h"
#include<QMenu>
#include<QAction>
#include<QMenuBar>
#include<QKeySequence>
#include<QMouseEvent>
#include<QPixmap>
#include<QDebug>
#include<QPushButton>
#include<QPainter>
#include<QList>
#include<QDialog>
#include<QInputDialog>
//2334
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    openAction=new QAction(QIcon(),"&open",this);
    openAction->setShortcuts(QKeySequence::Open);
    QMenu *file=menuBar()->addMenu("File");
    file->addAction(openAction);

    ps=new Pipsys(this);

    //
    elementsList=ps->elementsList;



    //Creat button group
    elementsBtns=new QButtonGroup(this);
    //Virtual button for that no element button is checked
    virtualBtn=new QPushButton(this);
    virtualBtn->setVisible(false);
    virtualBtn->setCheckable(true);
    elementsBtns->addButton(virtualBtn,0);
    //Creat action group
    elementActs=new QActionGroup(this);
    //Creat element actions
    addElementAction();

    //Allow unchecking all actions
    elementActs->setExclusionPolicy(QActionGroup::ExclusionPolicy::ExclusiveOptional);
    // Add element actions to toolbar
    ui->toolBar->addActions(elementActs->actions());
    //Creat pipe system object


}

MainWindow::~MainWindow()
{
    delete ui;
}
//[1]
//重写鼠标事件，用于在工作区生成元素按钮
void MainWindow::mousePressEvent(QMouseEvent *event){
    //Left button clicked in main window
    if(event->button()==Qt::LeftButton)
    {
        // Checked element action exists
        if(QAction *checkedAct= elementActs->checkedAction())
        {
            //Justify which element is activate
            QString &eleName=checkedAct->text();
            int eleID=elementsList.indexOf(eleName);
            if(eleID==1)//Pipe
            {
                ;
            }
            else
            {  // New button
                EleButton *button;
                button=new EleButton(ElementType(eleID),0,this);
                button->setGeometry(event->pos().x()-18,event->pos().y()-18,36,36);
                button->setCheckable(true);
                button->setStyleSheet(button_styleList.at(eleID));
//                connect(button,&QPushButton::d,this,&MainWindow::reservoirDialog);
                connect(button,&EleButton::doubleClicked,this,&MainWindow::showElementDialog);
                elementsBtns->addButton(button);
                ps->addElement(eleName);
                button->show();
            }
        }
        else{
            virtualBtn->toggle();
        }
        return;
    }
}

//[2]
//根据系统中的元素名称列表，生成工具栏动作
//初始元素按钮选中和非选中时的图标样式
void MainWindow::addElementAction()
{
    for (auto i = ps->elementsList.begin(); i != ps->elementsList.end(); ++i)
    {
        //New element actions according to elementsList
        QAction *elementAct=new QAction(QIcon(":/images/"+*i),*i,elementActs);
        elementAct->setCheckable(true);

        //New button styles
        QString button_style="QPushButton{\
                border-image:url(:/images/"+*i+")}"
                "QPushButton:checked {\
                border-image:url(:/images/"+*i+"_checked);}";
        button_styleList.append(button_style);
    }
};

//[3]
//分发元素对话框显示
void MainWindow::showElementDialog(int id, ElementType type)
{
    switch (type) {
    case ElementType::Reservoir:
        reservoirDialog(id);
        break;
    default:
        break;
    }
};
//[4]
//水库对话框
void MainWindow::reservoirDialog(int typeID)
{
      double waterLevel = QInputDialog::getDouble(this,
                         "Property",
                         "Water level" );
      ps->reservoirs[typeID].waterLevel=waterLevel;
};
