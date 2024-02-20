using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace PhysicShapes
{
    public class ThreeDBall
    {
        public Vector3 vel;
        public Vector3 pos;
        public float radius;
        PhysicScene s;

        public ThreeDBall(Vector3 vel, Vector3 pos, float radius)
        {
            this.vel = vel;
            this.pos = pos;
            this.radius = radius;
            s = PhysicScene.Instance;
        }

        public void simulate(float timestamp)
        {
            vel += timestamp * new Vector3(0, s.gravity, 0);
            vel *= .99f;
            pos += timestamp * vel;

            //����ǽ����ײ, todo ����λ�õķ�ʽ����Ӧ�ý���͸���־��룬�������µ�λ���ϣ���������ײ������ʱ�䣩
            if (pos.y < s.plane_lower_bound + radius)
            {
                pos.y = s.plane_lower_bound + radius; vel.y = -vel.y;
            }
            if (pos.x < s.x_lower_bound + radius)
            {
                pos.x = s.x_lower_bound + radius; vel.x = -vel.x;
            }
            if (pos.x > s.x_upper_bound - radius)
            {
                pos.x = s.x_upper_bound - radius; vel.x = -vel.x;
            }
            if (pos.z < s.z_lower_bound + radius)
            {
                pos.z = s.z_lower_bound + radius; vel.z = -vel.z;
            }
            if (pos.z > s.z_upper_bound - radius)
            {
                pos.z = s.z_upper_bound - radius; vel.z = -vel.z;
            }
        }
    }

    public class TwoDBall
    {
        public Vector2 vel;
        public float radius;
        public Vector2 pos;
        public float mass;
        public float restitution;  //������ײ�ָ�ϵ��

        public TwoDBall(Vector2 vel, float radius, Vector3 pos, float mass, float restitution)
        {
            this.vel = vel;
            this.radius = radius;
            this.pos = new Vector2(pos.x, pos.y);
            this.mass = mass;
            this.restitution = restitution;
        }
        public void Simulate(float gravity, float dt)
        {
            this.vel += new Vector2(0, -gravity) * dt;
            this.pos += this.vel * dt;
        }
    }
}