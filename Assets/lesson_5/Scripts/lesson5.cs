using PhysicShapes;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Bead : TwoDBall
{
    public Vector2 prevPos;
    public Bead(Vector2 vel, float radius, Vector3 pos, float mass, float restitution) : base(vel, radius, pos, mass, restitution) { }

    public void startStep(float dt, Vector2 gravity)
    {
        this.vel += dt * gravity;
        this.prevPos = this.pos;
        this.pos += dt * this.vel;
    }
    public void keepOnWire(Vector2 center, float radius)
    {
        var dir = this.pos - center;
        var len = dir.magnitude;
        if (len == 0.0)
            return;
        dir /= len;
        var lambda = radius - len;
        this.pos += lambda * dir;
    }
    public void endStep(float dt)
    {
        this.vel = (this.pos - this.prevPos) / dt;
    }
}

public class lesson5 : MonoBehaviour
{
    public Image img_ball1;
    public Image img_ball2;
    public Image img_ball3;
    Bead ball1;
    Bead ball2;
    Bead ball3;

    List<Bead> list_bead;

    public int numSteps = 1;
    void Start()
    {
        ball1 = new Bead(Vector2.zero, 75, img_ball1.transform.position, 75, 1);
        ball2 = new Bead(Vector2.zero, 50, img_ball2.transform.position, 50, 1);
        ball3 = new Bead(Vector2.zero, 50, img_ball3.transform.position, 50, 1);
        list_bead = new List<Bead>()
        {
            ball1, ball2, ball3
        };
    }

    public void BallBallCollsion(TwoDBall ball1, TwoDBall ball2)
    {
        var restitution = Mathf.Min(ball1.restitution, ball2.restitution);
        var dir = new Vector2();
        dir = ball2.pos - ball1.pos;
        var d = dir.magnitude;
        if (d == 0.0 || d > ball1.radius + ball2.radius) //去掉除0问题
            return;

        dir /= d;

        float corr = (ball1.radius + ball2.radius - d) / 2;
        ball1.pos -= corr * dir;
        ball2.pos += corr * dir;

        float v1 = Vector2.Dot(ball1.vel, dir);
        float v2 = Vector2.Dot(ball2.vel, dir);

        var m1 = ball1.mass;
        var m2 = ball2.mass;

        var newV1 = (m1 * v1 + m2 * v2 - m2 * (v1 - v2) * restitution) / (m1 + m2);
        var newV2 = (m1 * v1 + m2 * v2 - m1 * (v2 - v1) * restitution) / (m1 + m2);

        ball1.vel += (newV1 - v1) * dir; //碰撞方向速度归0所以去掉v1，在赋予新速度newV1。
        ball2.vel += (newV2 - v2) * dir;
    }

    private void FixedUpdate()
    {
        float dt = Time.fixedDeltaTime;

        var sdt = dt / numSteps;

        for (var step = 0; step < numSteps; step++)
        {
            for (var i = 0; i < list_bead.Count; i++)
                list_bead[i].startStep(sdt, Vector2.down * 9.8f * 200);

            for (var i = 0; i < list_bead.Count; i++)
            {
                list_bead[i].keepOnWire(new Vector2(Screen.width/2, Screen.height/2), 250f);
            }

            for (var i = 0; i < list_bead.Count; i++)
                list_bead[i].endStep(sdt);

            for (var i = 0; i < list_bead.Count; i++)
            {
                for (var j = 0; j < i; j++)
                {
                    BallBallCollsion(list_bead[i], list_bead[j]);
                }
            }
        }


        img_ball1.transform.position = ball1.pos;
        img_ball2.transform.position = ball2.pos;
        img_ball3.transform.position = ball3.pos;
    }
}
