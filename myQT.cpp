#include"myQT.h"
#include<QMessageBox>
EleButton::EleButton(ElementType type, int typeID,QWidget *parent)
    : QPushButton(parent)
{
    this->type=type;
    this->typeID=typeID;
}

EleButton::~EleButton()
{
}

void EleButton::mouseDoubleClickEvent(QMouseEvent * event)
{
    if (event->buttons() == Qt::LeftButton)
        emit doubleClicked(typeID,type);
}
