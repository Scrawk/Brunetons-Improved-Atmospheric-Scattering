using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace BrunetonsImprovedAtmosphere
{
    public class RotateCamera : MonoBehaviour
    {

        public float speed = 1.0f;

        private float viewDistanceMeters = 9.0f;

        private float viewZenithAngleRadians = 1.47f;

        private float viewAzimuthAngleRadians = -1.47f;

        private Vector3 lastMousePos;

        void Update()
        {
            HandleMouseWheelEvent();
            HandleMouseDragEvent();

            transform.LookAt(Vector3.zero);
        }

        private void HandleMouseWheelEvent()
        {
            float scroll = Input.GetAxis("Mouse ScrollWheel");

            if (scroll > 0.0f)
            {
                viewDistanceMeters *= 1.05f;
            }
            else if (scroll < 0.0f)
            {
                viewDistanceMeters /= 1.05f;
            }
        }

        private void HandleMouseDragEvent()
        {
            Vector3 delta = (lastMousePos - Input.mousePosition) * Time.deltaTime * speed;

            if (Input.GetMouseButton(0) && !Input.GetKey(KeyCode.LeftControl))
            {
                viewZenithAngleRadians += -delta.y;
                viewZenithAngleRadians = Mathf.Max(0.0f, Mathf.Min(Mathf.PI / 2.0f, viewZenithAngleRadians));
                viewAzimuthAngleRadians += delta.x;
            }

            float cos_z = Mathf.Cos(viewZenithAngleRadians);
            float sin_z = Mathf.Sin(viewZenithAngleRadians);
            float cos_a = Mathf.Cos(viewAzimuthAngleRadians);
            float sin_a = Mathf.Sin(viewAzimuthAngleRadians);

            Vector3 pos = new Vector3(sin_z * cos_a, cos_z, sin_z * sin_a) * viewDistanceMeters;

            transform.position = pos;

            lastMousePos = Input.mousePosition;
        }

    }
}
