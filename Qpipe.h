#ifndef QPIPE_H
#define QPIPE_H


#include <QGraphicsLineItem>

class Hydroment;

//! [0]
class QPipe : public QGraphicsLineItem
{
public:
    enum { Type = UserType + 4 };

    QPipe(Hydroment *startHydro, Hydroment *endHydro,
          QGraphicsItem *parent = nullptr);

    int type() const override { return Type; }
    QRectF boundingRect() const override;
    QPainterPath shape() const override;
    void setColor(const QColor &color) { myColor = color; }
    Hydroment *startItem() const { return myStartHydro; }
    Hydroment *endItem() const { return myEndHydro; }

    void updatePosition();

protected:
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
               QWidget *widget = nullptr) override;

private:
    Hydroment *myStartHydro;
    Hydroment *myEndHydro;
    QPolygonF QPipeHead;
    QColor myColor = Qt::black;
};
//! [0]


#endif // QPIPE_H
