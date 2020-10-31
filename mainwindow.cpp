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



    //Creat groups
    elementsBtns=new QButtonGroup(this);
    virtualBtn=new QPushButton(this);
    virtualBtn->setVisible(false);
    virtualBtn->setCheckable(true);
    elementsBtns->addButton(virtualBtn,0);
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
                QPushButton *button;
                button=new QPushButton(this);
                button->setGeometry(event->pos().x()-18,event->pos().y()-18,36,36);
                button->setCheckable(true);
                button->setStyleSheet(button_styleList.at(eleID));
                connect(button,&QPushButton::,this,&MainWindow::reservoirDialog);
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


void MainWindow::addElementAction()
{
    for (auto i = elementsList.begin(); i != elementsList.end(); ++i)
    {
        //New element actions
        QAction *elementAct=new QAction(QIcon(":/images/"+*i),*i,elementActs);
        elementAct->setCheckable(true);
//        //New button Groups
//        QButtonGroup* eleButtonGroup=new QButtonGroup(this);
//        eleBGList.append(eleButtonGroup);
        //New button styles
        QString button_style="QPushButton{\
                border-image:url(:/images/"+*i+")}"
                "QPushButton:checked {\
                border-image:url(:/images/"+*i+"_checked);}";
        button_styleList.append(button_style);
    }
};




void MainWindow::reservoirDialog(int typeID)
{
    QDialog* dialog=new QDialog(this);
    dialog->show();

};
