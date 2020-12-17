#ifndef MYQT_H
#define MYQT_H
#include<QPushButton>
#include<QMouseEvent>
#include<elements.h>
#include<QAction>
enum class ElementType
{
    Reservoir = 0,
    Pipe = 1,
    Node = 2,
    Branch = 3,
    Turbine = 4,
    BallValve = 5,
    EndBallValve = 6,
    Pump = 7,
    EndValve = 8,
    number = 9
};

class EleButton:public QPushButton
{
    Q_OBJECT
public:
    int typeID;
    ElementType type;
    EleButton(ElementType type, int typeID=-1,QWidget *parent=0);

    ~EleButton();
//    void remove
signals:
    void doubleClicked(int id, ElementType type) const;
protected:
    void mouseDoubleClickEvent(QMouseEvent *event); /*双击事件响应函数*/
};



class EleAction:public QAction
{
    Q_OBJECT
public:

    ElementType type;
    EleAction(ElementType type, int typeID=-1,QWidget *parent=0);
    ~EleAction();

};

#endif // MYQT_H
